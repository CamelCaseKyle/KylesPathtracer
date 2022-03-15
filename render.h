//build passthrough into renderer

#pragma once

#define DBG_RENDER
#define GLM_FORCE_SWIZZLE

#include "shader.h"

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtx/quaternion.hpp>

#include <opencv2\core.hpp>

//forward declare here
struct Content;
class Context;


static std::unordered_map<int, std::vector<int>> textureHelper = {
	{ CV_8UC1, { GL_R8, GL_RED, GL_UNSIGNED_BYTE } },
	{ CV_8UC3, { GL_RGB, GL_BGR, GL_UNSIGNED_BYTE } },
	{ CV_8UC4, { GL_RGBA, GL_BGRA, GL_UNSIGNED_BYTE } },
	{ CV_16UC1, { GL_R16, GL_RED, GL_UNSIGNED_SHORT } },
	{ CV_32FC1, { GL_R32F, GL_RED, GL_FLOAT } }
};


//creates a texture on the GPU
GLuint createTexture(GLuint mag, GLuint min, GLuint wrap, bool unbind = true);
//creates a framebuffer with depth buffer
GLuint createFrameBuffer(int width, int height, GLuint& color_attachment, GLuint& depth_attachment, bool unbind = true);
//creates a basic framebuffer
GLuint createFrameBuffer(int width, int height, GLuint& color_attachment, bool unbind = true);
//creates a buffer on the GPU
GLuint createBuffer();
//creates a basic vert array
GLuint createVertexArray(bool unbind = true);

//uploads a texture to the GPU
void uploadTexture(cv::Mat& texture, GLuint textureID, bool unbind = true);
//uploads a vertex array to the GPU
void uploadBuffer(GLuint type, const void* buffer, GLuint size, GLuint vboID, bool unbind = true);
//get texture from GPU
void downloadTexture(cv::Mat& texture, GLuint textureID, bool unbind = true);

//fullscreen quad, indexed
const float passthrough_verts[8] = {
	-1.0f, -1.0f,
	-1.0f,  1.0f,
	 1.0f, -1.0f,
	 1.0f,  1.0f,
};
const float passthrough_uvs[8] = {
	0.0f, 0.0f,
	0.0f, 1.0f,
	1.0f, 0.0f,
	1.0f, 1.0f,
};
const unsigned int passthrough_inds[4] = {
	0, 1, 2, 3
};

//this is a render context, todo: stream context
class Context {
private:
	//changing this alone wont do anything, need to use the method (same with aspect and FOV (except raytracer))
	glm::vec4 clearColor = glm::vec4(0.f, 0.f, 0.f, 1.f);
	//special passthrough shader
	Shader passthrough;

public:
	//index helpers
	const enum fboIDs	  { f_frame };
	const enum rentexIDs  { r_frame, r_component };
	const enum shaderIDs  { s_passthrough };
	const enum vaoIDs	  { a_quad };
	const enum vboIDs	  { v_quad, v_quad_uv, v_quad_ind };

	Shader *shaders[8];
	//this is all the renderer specific content
	GLuint vao[8], fbo[8], resID[8], textures[8], vbo[16];
	//handle for window
	GLFWwindow *window = nullptr;
	//why not
	std::string name;
	//projection matrix
	glm::mat4 projection = glm::mat4(1.f);
	//viewport for screen
	glm::vec4 viewp = glm::vec4(0.f, 0.f, 0.f, 0.f);
	//aspect ratio and fov
	float aspect = 0.f, fov = 0.f, globalTime = 0.f;
	//mode, global time
	unsigned int fnum = 0;
	//handle for this context
	int ID = -1;

	//init glfw window
	Context(int width, int height, float _fov, char *title);
	//init passthrough shader
	bool initPassthrough();
	
	//sets the clear color
	void ClearColor(glm::vec4& newCol);
	//shouldnt be called
	void Reshape(float _aspect, float _fov);

	//show on screen
	void display(const GLuint target);
};


//The base container for openGL content
struct Content {
	//render context
	Context* target = nullptr;
	//shader to use
	Shader* shade = nullptr;
	//Rendering attributes and data
	const float *verts = nullptr, *norms = nullptr, *uvs = nullptr, *attr1 = nullptr, *attr2 = nullptr;
	const unsigned int *inds = nullptr, *attr3 = nullptr;
	//content model matrix
	glm::mat4 model = glm::mat4(1.f), 
		//view matrix is from tracking usually
		view = glm::mat4(1.f), 
		//inverse view transform for transforming camera instead of objects
		view_ray = glm::mat4(1.f), 
		//normal matrix is transpose(inverse(MV))
		normalMatrix = glm::mat4(1.f);
	//opengl resources (frame buffers, transform buffers, render textures, textures, vertex arrays, vertex buffers)
	GLuint glFBO[8], glTBO[4], glRenTex[4], glTex[8], glVAO[8], glVBO[16];
	//yep
	int n_verts = 0, n_uvs = 0, n_inds = 0;
};
