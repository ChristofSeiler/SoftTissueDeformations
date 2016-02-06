/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	GPU implemenation of itkMRFRegularizationFilter.
 */

#ifndef _itkMRFRegularizationFilterGPU_txx
#define _itkMRFRegularizationFilterGPU_txx
#include "itkMRFRegularizationFilterGPU.h"

#include "itkImageRegionIterator.h"
#include "itkConstantBoundaryCondition.h"

#include <iostream>
#include <fstream>

namespace itk
{

template<class TInputImage, class TOutputImage, class TNLabelImage>
MRFRegularizationFilterGPU<TInputImage,TOutputImage,TNLabelImage>
::MRFRegularizationFilterGPU(void):
  m_NumberOfIterations(40), m_GPUType(0)
{

  if( (int)InputImageDimension != (int)OutputImageDimension )
    {
    OStringStream msg;
    msg << "Input image dimension: " << InputImageDimension << 
		" != output image dimension: " 
		<< OutputImageDimension; 
    throw ExceptionObject(__FILE__, __LINE__,msg.str().c_str(),ITK_LOCATION);
    }
}

template<class TInputImage, class TOutputImage, class TNLabelImage>
MRFRegularizationFilterGPU<TInputImage,TOutputImage,TNLabelImage>
::~MRFRegularizationFilterGPU(void)
{
}

template<class TInputImage, class TOutputImage, class TNLabelImage>
void
MRFRegularizationFilterGPU<TInputImage,TOutputImage,TNLabelImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent <<" MRF Regulariziation filter object " << std::endl;

  os << indent <<" Number of iterations: " << 
    m_NumberOfIterations << std::endl;

}// end PrintSelf

template<class TInputImage, class TOutputImage, class TNLabelImage>
void
MRFRegularizationFilterGPU<TInputImage,TOutputImage,TNLabelImage>
::GenerateData()
{
	// set input image to output image
	InputImageType::ConstPointer input = this->GetInput();
	OutputImageType::Pointer output = this->GetOutput();
	output->SetRegions(input->GetRequestedRegion());
	output->Allocate();
	ImageRegionConstIterator<InputImageType> in(input, input->GetRequestedRegion());
	ImageRegionIterator<OutputImageType> out(output, input->GetRequestedRegion());
	in.GoToBegin(); out.GoToBegin();
	for(in.GoToBegin(), out.GoToBegin();!in.IsAtEnd(); ++out, ++in)
		out.Set(in.Get());

	// yeah we are about to dive into the new world of GPU programming \\o o//
	// write to buffer and then to texture
	// start shading application again
	// do 40 times
	// then convert to itkImage and write to ouput pointer of this filter
	// done

	// ati
	rect_ati_rgba_32.name		= "TEXRECT - float_ATI - RGBA - 32";
    rect_ati_rgba_32.texTarget		= GL_TEXTURE_RECTANGLE_ARB;
    rect_ati_rgba_32.texInternalFormat	= GL_RGBA_FLOAT32_ATI;
    rect_ati_rgba_32.texFormat		= GL_RGBA;
	// nvidia
	rect_nv_rgba_32.name		= "TEXRECT - float_NV - RGBA - 32";
    rect_nv_rgba_32.texTarget		= GL_TEXTURE_RECTANGLE_ARB;
    rect_nv_rgba_32.texInternalFormat	= GL_FLOAT_RGBA32_NV;
    rect_nv_rgba_32.texFormat		= GL_RGBA;

	// struct actually being used
	struct_textureParameters textureParameters;
	if(m_GPUType == 0)
		textureParameters = rect_ati_rgba_32;
	else
		textureParameters = rect_nv_rgba_32;

    // declare texture size, the actual data will be a vector 
    // of size texSize*texSize*4
    // create test data and fill arbitrarily
	int texSize = output->GetRequestedRegion().GetSize()[0];
    float* data = (float*)malloc(4*texSize*texSize*sizeof(float));
	float* dataConst = (float*)malloc(4*texSize*texSize*sizeof(float));
	float* dataConst2 = (float*)malloc(4*texSize*texSize*sizeof(float));
    float* result = (float*)malloc(4*texSize*texSize*sizeof(float));
	unsigned int count = 0;
	ImageRegionConstIterator<LabelImageType> labelIter(m_Labels, m_Labels->GetRequestedRegion());
	ImageRegionConstIterator<VectorImageType> obsIter(m_ObservationImage, 
		m_ObservationImage->GetRequestedRegion());
	ImageRegionConstIterator<VectorImageType> confidenceIter(m_ConfidenceImage, 
		m_ConfidenceImage->GetRequestedRegion());
	for(out.GoToBegin(), labelIter.GoToBegin(), obsIter.GoToBegin(), confidenceIter.GoToBegin(); 
			!out.IsAtEnd(); ++out, ++labelIter, ++obsIter, ++confidenceIter) {
		OutputImageType::PixelType value = out.Get();
		OutputImageType::PixelType obsValue = obsIter.Get();
		OutputImageType::PixelType confidenceValue = confidenceIter.Get();
		data[count] = value[0];
		dataConst[count] = obsValue[0];
		dataConst2[count] = confidenceValue[0];
		++count;
		data[count] = value[1];
		dataConst[count] = obsValue[1];
		dataConst2[count] = confidenceValue[1];
		++count;
		data[count] = 0;
		dataConst[count] = labelIter.Get();
		dataConst2[count] = 0;
		++count;
		data[count] = 0;
		dataConst[count] = 0;
		dataConst2[count] = 0;
		++count;
	}

    // viewport transform for 1:1 pixel=texel=data mapping
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0,texSize,0.0,texSize);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0,0,texSize,texSize);
	checkGLErrors ("glViewport");

    // create FBO and bind it (that is, use offscreen render target)
    GLuint fb;
    glGenFramebuffersEXT(1,&fb); 
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fb);
	checkGLErrors ("glBindFramebufferEXT");

	// two textures identifiers referencing old and new
	GLuint texID[2];
	// ping pong management vars
	int writeTex = 0;
	int readTex = 1;
	GLenum attachmentpoints[] = { GL_COLOR_ATTACHMENT0_EXT, 
							  GL_COLOR_ATTACHMENT1_EXT 
							};
	// attach two textures to FBO
	glGenTextures (2, texID);
	checkGLErrors ("glGenTextures");

    // create texture 0
	createTexture(texID[writeTex], texSize, textureParameters);

	// create texture 1
	createTexture(texID[readTex], texSize, textureParameters);
	transferToTexture(data, texID[readTex], texSize, textureParameters);

	// create texture for read-only data
	GLuint readTexID;
	glGenTextures (1, &readTexID);
	createTexture(readTexID, texSize, textureParameters);
	transferToTexture(dataConst, readTexID, texSize, textureParameters);

	GLuint readTexID2;
	glGenTextures (1, &readTexID2);
	createTexture(readTexID2, texSize, textureParameters);
	transferToTexture(dataConst2, readTexID2, texSize, textureParameters);

	// attach two textures to FBO
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachmentpoints[writeTex], textureParameters.texTarget, 
		texID[writeTex], 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachmentpoints[readTex], textureParameters.texTarget, 
		texID[readTex], 0);

	// set parameter
	cgSetParameter2f(dimension, texSize, texSize);
	// enable label and observation texture (read-only, not changed during the iteration)
    cgGLSetTextureParameter(observation, readTexID);
    cgGLEnableTextureParameter(observation);
	cgGLSetTextureParameter(confidence, readTexID2);
    cgGLEnableTextureParameter(confidence);

	for (int i = 0; i < m_NumberOfIterations; i++) {
		// set render destination
		glDrawBuffer (attachmentpoints[writeTex]);
		// enable texture (read-only)
		cgGLSetTextureParameter(image, texID[readTex]);
		cgGLEnableTextureParameter(image);

		// and render multitextured viewport-sized quad
		// make quad filled to hit every pixel/texel
		glPolygonMode(GL_FRONT,GL_FILL);
		// and render quad
		glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0); 
		glVertex2f(0.0, 0.0);
		glTexCoord2f(texSize, 0.0); 
		glVertex2f(texSize, 0.0);
		glTexCoord2f(texSize, texSize); 
		glVertex2f(texSize, texSize);
		glTexCoord2f(0.0, texSize); 
		glVertex2f(0.0, texSize);
		glEnd();
		// swap role of the two textures (read-only source becomes 
		// write-only target and the other way round)
		int tmp = writeTex;
		writeTex = readTex;
		readTex = tmp;

		// debug
		//if(texSize != 1) break;
	}
	checkGLErrors ("glEnd");

	// and read back
    glReadBuffer(attachmentpoints[readTex]);
    glReadPixels(0, 0, texSize, texSize,textureParameters.texFormat,GL_FLOAT,result);
	
	/*
    // print out results
    printf("Data before roundtrip:\n");
    for (int i=0; i<texSize*texSize*4; i++)
        printf("%f\n",data[i]);
    printf("Data after roundtrip:\n");
    for (int i=0; i<texSize*texSize*4; i++)
        printf("%f\n",result[i]);
	*/

	count = 0;
	for(out.GoToBegin(); !out.IsAtEnd(); ++out) {
		OutputImageType::PixelType value = out.Get();
		value[0] = result[count];
		value[1] = result[++count];
		out.Set(value);
		++count;
		++count;
		++count;
	}

    // clean up
	glDeleteFramebuffersEXT(1,&fb);
    free(data);
    free(result);
    glDeleteTextures (2,texID);
	glDeleteTextures (1,&readTexID);
	glDeleteTextures (1,&readTexID2);

}// end GenerateData

// callback function
void cgErrorCallback(void) {
    CGerror lastError = cgGetError();
    if(lastError) {
        printf(cgGetErrorString(lastError));
        printf(cgGetLastListing(cgContext));
    }
}

/**
 * Checks for OpenGL errors.
 * Extremely useful debugging function: When developing, 
 * make sure to call this after almost every GL call.
 */
void checkGLErrors (const char *label) {
    GLenum errCode;
    const GLubyte *errStr;
    
    if ((errCode = glGetError()) != GL_NO_ERROR) {
	errStr = gluErrorString(errCode);
	printf("OpenGL ERROR: ");
	printf((char*)errStr);
	printf("(Label: ");
	printf(label);
	printf(")\n.");
    }
}

/**
 * Checks framebuffer status.
 * Copied directly out of the spec, modified to deliver a return value.
 */
bool checkFramebufferStatus() {
    GLenum status;
    status = (GLenum) glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    switch(status) {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
            return true;
        case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
	    printf("Framebuffer incomplete, incomplete attachment\n");
            return false;
        case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
	    printf("Unsupported framebuffer format\n");
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
	    printf("Framebuffer incomplete, missing attachment\n");
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
	    printf("Framebuffer incomplete, attached images must have same dimensions\n");
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
	    printf("Framebuffer incomplete, attached images must have same format\n");
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
	    printf("Framebuffer incomplete, missing draw buffer\n");
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
	    printf("Framebuffer incomplete, missing read buffer\n");
            return false;
    }
    return false;
}

void initGPU() {
	// set up glut to get valid GL context and 
    // get extension entry points
	int argc = 1;
	char *argv = "bla.exe";
    glutInit (&argc, &argv);
    glutWindowHandle = glutCreateWindow("TEST1");
	int err = glewInit();
    // sanity check
    if (GLEW_OK != err) {
		printf((char*)glewGetErrorString(err));
    }  

	// set up Cg
	cgContext = cgCreateContext();
	// register the error callback once the context has been created
	cgSetErrorCallback(cgErrorCallback);
	fragmentProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
	cgGLSetOptimalOptions(fragmentProfile);
	// create fragment program
	fragmentProgram = cgCreateProgramFromFile(cgContext, CG_SOURCE, "./src/shader.cg",
		fragmentProfile, "energy_function", NULL);
	// load program
	cgGLLoadProgram (fragmentProgram);
	image = cgGetNamedParameter (fragmentProgram,"image");
	dimension = cgGetNamedParameter (fragmentProgram,"dimension");
	observation = cgGetNamedParameter (fragmentProgram,"observation");
	confidence = cgGetNamedParameter (fragmentProgram,"confidence");

	// enable fragment profile
	cgGLEnableProfile(fragmentProfile);
	// bind program
	cgGLBindProgram(fragmentProgram);
}

void createTexture(GLuint texID, int texSize, struct_textureParameters& textureParameters) {
	glBindTexture(textureParameters.texTarget,texID);
	checkGLErrors ("glBindTexture");
    // set texture parameters
    glTexParameteri(textureParameters.texTarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(textureParameters.texTarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(textureParameters.texTarget, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(textureParameters.texTarget, GL_TEXTURE_WRAP_T, GL_CLAMP);
	checkGLErrors ("glTexParameteri");
    // define texture with floating point format
    glTexImage2D(textureParameters.texTarget,0,textureParameters.texInternalFormat,texSize,texSize,0,
		textureParameters.texFormat,GL_FLOAT,0);
	checkGLErrors ("glTexImage2D");
}

void transferToTexture(float* data, GLuint texID, int texSize, 
					   struct_textureParameters& textureParameters) {
    // transfer data to texture
	if(textureParameters.texInternalFormat == GL_RGBA_FLOAT32_ATI) {
		// version (b): HW-accelerated on ATI 
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
			textureParameters.texTarget, texID, 0);
		glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
		checkGLErrors ("glDrawBuffer");
		glRasterPos2i(0,0);
		glDrawPixels(texSize,texSize,textureParameters.texFormat,GL_FLOAT,data);
		checkGLErrors ("glDrawPixels");
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
			textureParameters.texTarget, 0, 0);
	}
	else {
		// version (a): HW-accelerated on NVIDIA 
		glBindTexture(textureParameters.texTarget,texID);
		glTexSubImage2D(textureParameters.texTarget,0,0,0,texSize,texSize,
			textureParameters.texFormat,GL_FLOAT,data);
		checkGLErrors ("glTexSubImage2D");
	}
}

void cleanUpGPU() {
	cgDestroyProgram(fragmentProgram);
    cgDestroyContext(cgContext);
	// will crash in debug mode on NVIDIA
    //glutDestroyWindow (glutWindowHandle);
}

} // end namespace itk

#endif