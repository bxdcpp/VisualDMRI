#ifndef TRKFILEIO_H
#define TRKFILEIO_H
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "defs.h"
using namespace std;

static const int TRK_HEADER_SIZE = 1000;    ///< header size is 1000
#pragma pack(push)                          ///< push current alignment to stack
#pragma pack(1)                             ///< set alignment to 1 byte boundary
/**
 * @brief The TrkFileHeader class file header defined at http://www.trackvis.org/docs/?subsect=fileformat
 */
class TrkFileHeader
{
public:
    TrkFileHeader()
    {
        memset(this, 0, TRK_HEADER_SIZE);
        strcpy_s(id_string,"TRACK");
        strcpy_s(voxel_order,"LAS");
        version = 2;
        hdr_size = 1000;
    }
    char id_string[6];                          /// char	6	ID string for track file. The first 5 characters must be "TRACK".
    int16_t dim[3];                             /// short int	6	Dimension of the image volume.
    float voxel_size[3] ;                       /// float	12	Voxel size of the image volume.
    float origin[3];                            ///	float	12	Origin of the image volume. This field is not yet being used by TrackVis. That means the origin is always (0, 0, 0).
    int16_t n_scalars;                          ///	short int	2	Number of scalars saved at each track point (besides x, y and z coordinates).
    char scalar_name[10][20];                   ///	char	200	Name of each scalar. Can not be longer than 20 characters each. Can only store up to 10 names.
    int16_t n_properties;                       ///	short int	2	Number of properties saved at each track.
    char property_name[10][20];                 ///	char	200	Name of each property. Can not be longer than 20 characters each. Can only store up to 10 names.
    float vox_to_ras[4][4];                     ///	float	64	4x4 matrix for voxel to RAS (crs to xyz) transformation. If vox_to_ras[3][3] is 0, it means the matrix is not recorded. This field is added from version 2.
    char reserved[444];                         ///	char	444	Reserved space for future version.
    char voxel_order[4];                        ///	char	4	Storing order of the original image data. Explained here.
    char pad2[4];                               ///	char	4	Paddings.
    float image_orientation_patient[6];         ///	float	24	Image orientation of the original image. As defined in the DICOM header.
    char pad1[2];                               ///	char	2	Paddings.
    unsigned char invert_x;                     ///	unsigned char	1	Inversion/rotation flags used to generate this track file. For internal use only.
    unsigned char invert_y;                     ///	unsigned char	1	As above.
    unsigned char invert_z;                     ///	unsigned char	1	As above.
    unsigned char swap_xy;                      ///	unsigned char	1	As above.
    unsigned char swap_yz;                      ///	unsigned char	1	As above.
    unsigned char swap_zx;                      ///	unsigned char	1	As above.
    int32_t n_count;                            ///	int	4	Number of tracks stored in this track file. 0 means the number was NOT stored.
    int32_t version;                            ///	int	4	Version number. Current version is 2.
    int32_t hdr_size;                           ///	int	4	Size of the header. Used to determine byte swap. Should be 1000.
};
#pragma pack(pop)                               /// restore original alignment from stack


/**
 * @brief The TrkInfo struct Track offset, bytes and point numbers
 */
struct TrkInfo
{
    streamoff trackOffset;  ///< offset in file
    int16_t numPntsInTrk;   ///< points in this track
    int16_t lengthInByte;   ///< length in bytes for all points in this track
};

/**
 * @brief The TrkFileReader class track (*.trk) file writer
 */
class TrkFileReader
{
public:
    TrkFileReader()
		{
			xInit();
		}

    TrkFileReader(string& strFilepath)
		{
			xInit();
			m_cFilepath = strFilepath;
		}

    ~TrkFileReader()
		{
			close();
		}

    /**
     * @brief open Open the track file
     * @return
     */
    bool open()
		{
			/// open trk file
			m_cFile.open(m_cFilepath.c_str(), ios::in | ios::binary);
			if (!m_cFile.is_open())
			{
				std::cerr << "fail to open file" << std::endl;
				return false;
			}
			m_cFile.seekg(0);
			assert(sizeof(m_cHeader) == TRK_HEADER_SIZE);         ///< make sure header is 1000 bytes
			m_cFile.read((char*)(&m_cHeader), TRK_HEADER_SIZE);   ///< read header
			assert(m_cFile.gcount() == TRK_HEADER_SIZE);          ///< read check

			if (m_cHeader.n_properties != 0)
				std::cerr << "n_properties is not 0:" << m_cHeader.n_properties << std::endl;
			if (m_cHeader.n_scalars != 0)
				std::cerr << "n_scalars is not 0:" << m_cHeader.n_scalars << std::endl;

			/// build the random access table (track number to file offset map)
			static TrkInfo cTrkInfo;
			m_iTrkPos = 0;
			m_cRandomAccessMap.clear();
			int32_t iTrkSizeInByte = 0;
			while (!m_cFile.eof())
			{
				/// read track info
				m_cFile.read((char*)(&m_iPntPosMax), sizeof(int32_t));
				iTrkSizeInByte = ((3 + m_cHeader.n_scalars)*m_iPntPosMax + m_cHeader.n_properties) * sizeof(float);
				/// save track position and length in byte
				cTrkInfo.numPntsInTrk = m_iPntPosMax;
				cTrkInfo.trackOffset = m_cFile.tellg();
				cTrkInfo.lengthInByte = iTrkSizeInByte;
				m_cRandomAccessMap[m_iTrkPos] = cTrkInfo;
				/// move to next track
				m_cFile.seekg(iTrkSizeInByte, ios::cur);
				m_cFile.peek();
				++m_iTrkPos;
			}
			xResetPos();
			return true;
		}

    /**
     * @brief close Close the track file
     */
    void close()
		{
			xInit();
			m_cFile.close();
		}

    /**
     * @brief readTrack Read one track
     * @param iTrkIdx the track index (start from 0)
     * @param points
     * @return
     */
    bool readTrack(size_t iTrkIdx, vector<float>& points)
		{
			map<int32_t, TrkInfo>::iterator it = m_cRandomAccessMap.find(static_cast<int32_t>(iTrkIdx));
			if (it == m_cRandomAccessMap.end())
				return false;
			TrkInfo& cTrkInfo = it->second;
			int iTotalPoints = cTrkInfo.numPntsInTrk;
			static float afPoint[3];
			m_cFile.seekg(cTrkInfo.trackOffset);    /// seek to track
			points.clear();
			for (int i = 0; i < iTotalPoints; i++)
			{
				m_cFile.read((char*)(afPoint), sizeof(float) * 3);    /// read three component
				m_cFile.seekg(sizeof(float)*m_cHeader.n_scalars, ios::cur); /// ignore scalars
				for (int j = 0; j < 3; ++j)
					points.push_back(afPoint[j]);
			}
			return true;
		}

    /**
     * @brief readPoint Read one point in the track
     * @param iTrkIdx the track index (start from 0)
     * @param iPntIdx the point index (start from 0)
     * @param point
     * @return
     */
    bool readPoint(int iTrkIdx, int iPntIdx, vector<float>& point)
		{
			map<int32_t, TrkInfo>::iterator it = m_cRandomAccessMap.find(iTrkIdx);
			if (it == m_cRandomAccessMap.end())
				return false;
			TrkInfo& rcTrk = it->second;
			streamoff iOffset = rcTrk.trackOffset + (3 + m_cHeader.n_scalars)*iPntIdx * sizeof(float);
			m_cFile.seekg(iOffset);

			float fCoord;
			point.clear();
			for (int i = 0; i < 3; i++)
			{
				m_cFile.read((char*)(&fCoord), sizeof(float));
				point.push_back(fCoord);
			}
			return true;
		}

    /**
     * @brief getTotalTrkNum Get total number of tracks in this file
     * @return total number of tracks in this file
     */
    size_t getTotalTrkNum()
		{
			return m_cRandomAccessMap.size();
		}

    /**
     * @brief getPointNumInTrk Get total number of points in the track
     * @param iIdx track index (start from 0)
     * @return total number of the points in this track
     */
    size_t getPointNumInTrk(int iIdx)
		{
			map<int32_t, TrkInfo>::iterator it = m_cRandomAccessMap.find(iIdx);
			if (it == m_cRandomAccessMap.end())
				return -1;
			return it->second.numPntsInTrk;
		}

    /**
     * @brief checkFile For internal use, check is the file is valid
     */
    void checkFile()
		{
			std::cerr << "Start Checking..." << std::endl;
			m_cFile.seekg(0);

			assert(sizeof(m_cHeader) == TRK_HEADER_SIZE);      ///< make sure header is 1000 bytes
			m_cFile.read((char*)(&m_cHeader), TRK_HEADER_SIZE);  ///< read header

			std::cerr << "n_count:" << m_cHeader.n_count << "\nn_scalars:" << m_cHeader.n_scalars << "\nn_properties:" << m_cHeader.n_properties << std::endl;

			m_iPntPos = 0;
			m_iTrkPos = 0;
			m_cFile.seekg(TRK_HEADER_SIZE);

			float points[3];

			m_iTrkPos = 0;
			while (!m_cFile.eof())
			{

				m_cFile.read((char*)(&m_iPntPosMax), sizeof(int32_t));
				std::cerr << "Track" << m_iTrkPos << ":" << m_iPntPosMax << "points" << std::endl;

				for (m_iPntPos = 0; m_iPntPos < m_iPntPosMax; ++m_iPntPos)
				{
					m_cFile.read((char*)points, 3 * sizeof(float));
					streamsize lReadSize = m_cFile.gcount();
					if (lReadSize != 3 * sizeof(float))
					{
						std::cerr << "Reading Fail for #:" << m_iPntPos << "Point" << std::endl;
						std::cerr << "Read Size is:" << lReadSize << std::endl;
						assert(lReadSize == 3 * sizeof(float));
					}
					m_cFile.seekg(m_cHeader.n_scalars * sizeof(float), ios::cur);
				}
				m_cFile.seekg(m_cHeader.n_properties * sizeof(float), ios::cur);
				m_cFile.peek();
				++m_iTrkPos;
			}
			std::cerr << "Total:" << m_iTrkPos << std::endl;
			std::cerr << "Checking Finished..." << std::endl;

			xResetPos();
		}

protected:
    void xInit()
		{
			m_cFilepath.clear();
			m_iPntPos = -1;
			m_iTrkPos = -1;
			m_iPntPosMax = -1;
		}

    void xResetPos()
		{
			m_iPntPos = -1;
			m_iTrkPos = -1;
			m_iPntPosMax = -1;
		}

    ADD_CLASS_FIELD_PRIVATE(fstream , cFile)                        ///< track file stream
    ADD_CLASS_FIELD(string, cFilepath, getFilepath, setFilepath)    ///< track file path
    ADD_CLASS_FIELD_NOSETTER(TrkFileHeader, cHeader, getHeader)     ///< track file header
    ADD_CLASS_FIELD(int32_t, iTrkPos, getTrkPos, setTrkPos)         ///< current track index
    ADD_CLASS_FIELD(int32_t, iPntPos, getPntPos, setPntPos)         ///< current point index
    ADD_CLASS_FIELD_NOSETTER(int32_t, iPntPosMax, getPntPosMax)     ///< total point number in current track
    ADD_CLASS_FIELD_NOSETTER(CONCATE(map<int32_t,TrkInfo>), cRandomAccessMap, getRandomAccessMap)    /// track (index to offset) map for random access support
};


/**
 * @brief The TrkFileWriter class track (*.trk) file reader
 */
class TrkFileWriter
{
public:
    TrkFileWriter() {	}
    TrkFileWriter(string &strFilepath)
		{
			m_cFilepath = strFilepath;
		}

    ~TrkFileWriter()
		{
			close();
		}

    /**
     * @brief create Create a new empty trk file for write. **MUST BE EMPTY**
     * @return
     */
    bool create()
		{
			/// open trk file
			m_cFile.open(m_cFilepath.c_str(), ios::out | ios::binary);
			if (!m_cFile.is_open())
			{
				std::cerr << "fail to open file" << std::endl;
				return false;
			}
			xWriteHeader();
			return true;
		}

    /**
     * @brief copyHeader Copy a header form the other file
     * @param other
     */
    void copyHeader(const TrkFileHeader& other)
		{
			memcpy(&m_cHeader, &other, TRK_HEADER_SIZE);
			m_cHeader.n_count = 0;
			m_cHeader.n_properties = 0; /// TODO for now it simply ignore n_properties & n_scalars
			m_cHeader.n_scalars = 0;
		}

    /**
     * @brief appendTrack Append a track to the end of the file
     * @param points
     * @return
     */
    bool appendTrack(vector<float>& points)
		{
			/// check if 3d points
			if (points.size() % 3 != 0)
			{
				std::cerr << "point number cannot divided by 3" << std::endl;
				return false;
			}
			/// append
			m_cFile.seekp(0, ios::end);

			/// write total point number in this track
			int32_t iTotalPoint = static_cast<int32_t>(points.size() / 3);
			m_cFile.write((char*)&iTotalPoint, sizeof(int32_t));

			/// write the points

			m_cFile.write((char*)points.data(), points.size() * sizeof(float));

			/// fill extra scalars with 0
			float fZero = 0;
			for (int i = 0; i < m_cHeader.n_scalars; i++)
				m_cFile.write((char*)&fZero, sizeof(float));


			/// fill extra properties with 0
			for (int i = 0; i < m_cHeader.n_properties; i++)
				m_cFile.write((char*)&fZero, sizeof(float));

			return true;
		}

    /**
     * @brief save Save the file immediately
     */
    void save()
		{
			xWriteHeader();
			return;
		}

    /**
     * @brief close Close the file
     */
    void close()
		{
			xWriteHeader();
			m_cFile.close();
		}

protected:
    void xWriteHeader()
		{
			m_cFile.seekp(0);
			assert(sizeof(m_cHeader) == TRK_HEADER_SIZE);         ///< make sure header is 1000 bytes
			m_cFile.write((char*)(&m_cHeader), TRK_HEADER_SIZE);   ///< write header
		}

    ADD_CLASS_FIELD_PRIVATE(fstream , cFile)    ///< track file stream
    ADD_CLASS_FIELD_NOSETTER(TrkFileHeader, cHeader, getHeader)     ///< track file header
    ADD_CLASS_FIELD(string, cFilepath, getFilepath, setFilepath)    ///< track file path
};

#endif // TRKFILEIO_H
