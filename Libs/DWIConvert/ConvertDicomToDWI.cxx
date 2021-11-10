 /*=========================================================================
 Computing (NAMIC), funded by the National Institutes of Health
 through the NIH Roadmap for Medical Research, Grant U54 EB005149.

 See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

  ***
  This program converts Diffusion weighted MR images in Dicom format into
  NRRD format.

 Assumptions:

 1) Uses left-posterior-superior (Dicom default) as default space for philips and siemens.
 This is the default space for NRRD header.
 2) For GE data, Dicom data are arranged in volume interleaving order.
 3) For Siemens data, images are arranged in mosaic form.
 4) For oblique collected Philips data, the measurement frame for the
 gradient directions is the same as the ImageOrientationPatient

 Reference materials:
 DICOM Data Dictionary: http://medical.nema.org/Dicom/2011/11_06pu.pdf
 =========================================================================*/

#include "DWIConverter.h"

 /** the real computation goes on in DWIConverter classes, of which
  * there is one for each manufacturer we encounter.
  */
#include "DWIConverterFactory.h"
#include "FSLDWIConverter.h"
#include "NRRDDWIConverter.h"


#include "DWIConvertLib.h"

int
main(int argc, char* argv[])
{
    // const std::string version = commandLine.getVersion();
    // BRAINSRegisterAlternateIO();
    const std::string conversionMode = "DicomToNrrd";
    const std::string fslNIFTIFile = "";
    std::string inputDicomDirectory = "C:/Users/Bxd/Desktop/testhead/5-ep2d_diff_mddw_20_p2";
    std::string inputVolume = "";
    std::string inputBValues;
    std::string inputBVectors;
    const std::string gradientVectorFile;
    double smallGradientThreshold = 0.2;
    bool fMRIOutput = false;
    bool allowLossyConversion = false;
    bool useIdentityMeaseurementFrame = false;

    //if manufacturer is siemens useBMatrixGradientDirections is true
    bool useBMatrixGradientDirections = true;
    std::string outputNiftiFile = "";
    std::string outputVolume = "C:/Users/Bxd/Desktop/testhead/dwi.nrrd";
    const std::string outputDirectory = "";
    const std::string outputBValues = "";
    const std::string outputBVectors = "";



    std::cout << "======= DWI Convert Tool Program=========" << std::endl;
    DWIConvert dWIConvert;

    if (conversionMode == "FSLToNrrd" && (!fslNIFTIFile.empty()))
    {
        inputVolume = fslNIFTIFile;
    }

    const std::string inputPathName = inputVolume.empty() ? inputDicomDirectory : inputVolume;

    dWIConvert.SetInputFileName(inputPathName);
    dWIConvert.setInputBValues(inputBValues);
    dWIConvert.setInputBVectors(inputBVectors);
    dWIConvert.setGradientVectorFile(gradientVectorFile);
    dWIConvert.setSmallGradientThreshold(smallGradientThreshold);

    dWIConvert.setfMRIOutput(fMRIOutput);
    dWIConvert.setAllowLossyConversion(allowLossyConversion);
    dWIConvert.setUseIdentityMeasurementFrame(useIdentityMeaseurementFrame);
    dWIConvert.setUseBMatrixGradientDirections(useBMatrixGradientDirections);

    if (!outputNiftiFile.empty())
    {
        outputVolume = outputNiftiFile;
    }

    dWIConvert.SetOutputFileName(outputVolume);
    dWIConvert.setOutputDirectory(outputDirectory);
    dWIConvert.setOutputBValues(outputBValues);
    dWIConvert.setOutputBVectors(outputBVectors);

    int result = dWIConvert.read();
    if (EXIT_SUCCESS == result)
    {
        return dWIConvert.write(outputVolume);
    }
    else
        return result;
}
