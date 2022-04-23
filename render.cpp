#include "render.h"

using namespace std;
using namespace glm;

// creates a texture on the GPU
GLuint createTexture(GLuint mag, GLuint min, GLuint wrap, bool unbind) {
	GLuint newID;
	glGenTextures(1, &newID);
	glBindTexture(GL_TEXTURE_2D, newID);
	// parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap);
	// if should leave bound
	if (unbind) glBindTexture(GL_TEXTURE_2D, 0);
	return newID;
}

// creates a framebuffer with depth buffer
GLuint createFrameBuffer(int width, int height, GLuint& color_attachment, GLuint& depth_attachment, bool unbind) {
	// path trace render texture
	color_attachment = createTexture(GL_NEAREST, GL_NEAREST, GL_CLAMP, false);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, 0);
	// fbo depth attachment
	depth_attachment = createTexture(GL_NEAREST, GL_NEAREST, GL_CLAMP, false);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	// The frame buffer handle
	GLuint framebuffer;
	glGenFramebuffers(1, &framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	// attach vertex render textures
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_attachment, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_attachment, 0);
	// tell it what to do with the buffers
	GLenum drawBuffers[2] = { GL_COLOR_ATTACHMENT0, GL_DEPTH_ATTACHMENT };
	// "2" is the size of DrawBuffers
	glDrawBuffers(2, drawBuffers);
	if (unbind) glBindFramebuffer(GL_FRAMEBUFFER, 0);
	return framebuffer;
}

// creates a basic framebuffer
GLuint createFrameBuffer(int width, int height, GLuint& color_attachment, bool unbind) {
	// path trace render texture
	color_attachment = createTexture(GL_NEAREST, GL_NEAREST, GL_CLAMP, false);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, 0);
	// The frame buffer handle
	GLuint framebuffer;
	// The raytrace frame buffer
	glGenFramebuffers(1, &framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	// attach vertex render textures
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_attachment, 0);
	// tell it what to do with the buffers
	GLenum drawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, drawBuffers);
	if (unbind) glBindFramebuffer(GL_FRAMEBUFFER, 0);
	return framebuffer;
}

// creates a buffer on the GPU
GLuint createBuffer() {
	GLuint newBuff;
	glGenBuffers(1, &newBuff);
	return newBuff;
}

// creates a basic vert array
GLuint createVertexArray(bool unbind) {
	GLuint newVAO;
	glGenVertexArrays(1, &newVAO);
	if (!unbind) glBindVertexArray(newVAO);
	return newVAO;
}

// uploads a texture to the GPU
//void uploadTexture(cv::Mat& texture, GLuint textureID, bool unbind) {
//	vector<int>& info = textureHelper[texture.type()];
//	// bind texture in question
//	glBindTexture(GL_TEXTURE_2D, textureID);
//	// upload
//	glTexImage2D(GL_TEXTURE_2D, 0, info[0], texture.cols, texture.rows, 0, info[1], info[2], texture.data);
//	if (unbind) glBindTexture(GL_TEXTURE_2D, 0);
//}

// uploads a vertex array to the GPU
void uploadBuffer(GLuint type, const void* buffer, GLuint size, GLuint vboID, bool unbind) {
	glBindBuffer(type, vboID);
	glBufferData(type, size, buffer, GL_STATIC_DRAW);
	if (unbind) glBindBuffer(type, 0);
}

// get texture from GPU
//void downloadTexture(cv::Mat& texture, GLuint textureID, bool unbind) {
//	vector<int>& info = textureHelper[texture.type()];
//	// bind texture in question
//	glActiveTexture(GL_TEXTURE0);
//	glBindTexture(GL_TEXTURE_2D, textureID);
//	// download
//	glGetTexImage(GL_TEXTURE_2D, 0, info[1], info[2], texture.data);
//	if (unbind) glBindTexture(GL_TEXTURE_2D, 0);
//}

// initialize GLFW, OpenGL
Context::Context(int width, int height, float _fov, char *title) {
	name = title;

	// Initialise GLFW
	if (!glfwInit()) {
		printf("Failed to initialize GLFW\n");
		name = nullptr;
		return;
	}

	//glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	glfwWindowHint(GLFW_FLOATING, GL_TRUE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	
	window = glfwCreateWindow(width, height, title, NULL, NULL);

	if (window == NULL) {
		printf("Failed to open GLFW window\n");
		glfwTerminate();
		name = nullptr;
		return;
	}

	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		printf("Failed to initialize GLEW\n");
		name = nullptr;
		return;
	}

	// ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	// section on screen
	viewp = glm::vec4(0, 0, width, height);
	// ratio
	aspect = float(width) / float(height);
	// fov
	fov = _fov;
	// Projection matrix : fov and aspect from camera, display range : 50cm - 4m
	projection = glm::perspective(fov, aspect, 0.1f, 1000.0f);
	// set later
	clearColor = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
	// all teh matricies
	glViewport(0, 0, int(viewp[2]), int(viewp[3]));
	// GL environment
	glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);

	// proper occlusion
	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LESS);
	//glClearDepth(1.0);
}

// initialize passthrough texture renderer
bool Context::initPassthrough() {
	passthrough = Shader("passthrough.vert", "passthrough.frag");
	passthrough.addUniforms({ "iResolution", "iFrame", "iChannel0", "iChannel1", "iChannel2"});
	shaders[s_passthrough] = &passthrough;

	vao[a_quad] = createVertexArray(false);

	// vert buffer
	vbo[v_quad] = createBuffer();
	uploadBuffer(GL_ARRAY_BUFFER, &passthrough_verts[0], 4 * 2 * sizeof(float), vbo[v_quad], false);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
	// uv buffer
	vbo[v_quad_uv] = createBuffer();
	uploadBuffer(GL_ARRAY_BUFFER, &passthrough_uvs[0], 4 * 2 * sizeof(float), vbo[v_quad_uv], false);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
	// indices
	vbo[v_quad_ind] = createBuffer();
	uploadBuffer(GL_ELEMENT_ARRAY_BUFFER, &passthrough_inds[0], 4 * sizeof(unsigned int), vbo[v_quad_ind], false);

	glBindVertexArray(0);

	// The pathtrace frame buffer
	fbo[f_frame] = createFrameBuffer(int(viewp[2]), int(viewp[3]), resID[r_frame], resID[r_component], false);

	// best check ever
	GLuint status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT ||
		status == GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT ||
		status == GL_FRAMEBUFFER_UNSUPPORTED)
	{
		printf("Failed to generate context %i\n", status);
		return false;
	}

	// unbind frame buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	return true;

}

// set the clear color
void Context::ClearColor(glm::vec4& newCol) {
	// keep track
	clearColor = newCol;
	// send to GL environment
	glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
}

// set window size
void Context::Reshape(float _aspect, float _fov) {
	// ratio
	aspect = _aspect;
	// view feild
	fov = _fov;
	// recalc projection
	projection = glm::perspective(fov, aspect, 50.0f, 400.0f);
}

// display everything on screen
void Context::display(const GLuint gBuff, const GLuint dBuff, const GLuint sBuff) {
	// render to screen
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// use our shader?
	if (shaders[s_passthrough] == nullptr) {
		printf("no passthrough shader\n");
		return;
	}
	Shader &shade = *shaders[s_passthrough];
	glUseProgram(shade.getID());

	// use vao
	glBindVertexArray(vao[a_quad]);

	// some uniforms...
	glUniform2f(shade.uniform("iResolution"), viewp[2], viewp[3]);
	glUniform1i(shade.uniform("iFrame"), fnum);
	glUniform1f(shade.uniform("iTime"), globalTime);
	//glUniform3f(shade.uniform("loc"), 0.0f, 0.0f, 0.0f);
	//glUniform3f(shade.uniform("vel"), 0.0f, 0.0f, 0.0f);
	//glUniform3f(shade.uniform("iMouse"), 0.0f, 0.0f, 0.0f);
	//glUniform3f(shade.uniform("orient"), 0.0f, 0.0f, 0.0f);
	
	// bind textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, gBuff);
	glUniform1i(shade.uniform("iChannel0"), 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, dBuff);
	glUniform1i(shade.uniform("iChannel1"), 1);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, sBuff);
	glUniform1i(shade.uniform("iChannel2"), 2);
	//glActiveTexture(GL_TEXTURE3);
	//glBindTexture(GL_TEXTURE_2D, resID[r_ray]);
	//glUniform1i(shade.uniform("iChannel3"), 3);

	// draw the indexed triangles!
	glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_INT, 0);

	//glDisable(GL_TEXTURE0);
	//glDisable(GL_TEXTURE1);
	//glDisable(GL_TEXTURE2);
	//glDisable(GL_TEXTURE3);
	//glDisable(gBuff);
	//glDisable(target);
}
