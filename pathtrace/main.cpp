#define GLM_SWIZZLE

#include <time.h>
#include <string>
#include <opencv2\core.hpp>
#include "render.h"

using namespace glm;

//index helpers
const enum fboIDs	  { f_frame };
const enum rentexIDs  { r_frame, r_component };
const enum shaderIDs  { s_pathtrace };
const enum vaoIDs	  { a_quad };
const enum vboIDs	  { v_quad, v_quad_uv, v_quad_ind };

const float quad_verts[8] = {
	-1.f, -1.f,
	-1.f, 1.f,
	1.f, -1.f,
	1.f, 1.f,
};
const float quad_uvs[8] = {
	-1.f, -1.f,
	-1.f, 1.f,
	1.f, -1.f,
	1.f, 1.f,
};
const unsigned int quad_inds[6] = {
	0, 1, 2,
	2, 1, 3,
};

const float ACCEL = 0.01f,
			LOOK = 0.005f,
			zfar = 1000.f,
			eps = .00001f, ieps = .99999f,
			sml = .001f, isml = .999f,
			sc45 = .7071067f,
			pi_5 = 1.5707963f, pi2 = 6.2831853f,
			pi = 3.1415926f, pi_rcp = .3183098f;

vec3 loc, vel, orient;
double mouseX, mouseY, mouseP, imx, imy;

inline float sign(float val) { return float((0.f < val) - (val < 0.f)); }

vec3 rotateY(const vec3& p, const vec2& angle) {
	vec2 c = cos(angle), s = sin(angle); vec3 o = p;
	//o.yz = o.yz * mat2(c.x, s.x, -s.x, c.x);
	o.xz = o.xz * mat2(c.y, s.y, -s.y, c.y);
	return o;
}

//init content
void initGeometry(Content& cont) {
	cont.model = mat4(1.0);
	cont.view = mat4(1.0);
	cont.view_ray = mat4(1.0);
	cont.normalMatrix = mat4(1.0);
	cont.verts = quad_verts;
	cont.uvs = quad_uvs;
	cont.inds = quad_inds;
	cont.n_verts = 4 * 2;
	cont.n_uvs = 4 * 2;
	cont.n_inds = 6;

	cont.glVAO[a_quad] = createVertexArray(false);

	//vert buffer
	cont.glVBO[v_quad] = createBuffer();
	uploadBuffer(GL_ARRAY_BUFFER, &cont.verts[0], cont.n_verts * sizeof(float), cont.glVBO[v_quad], false);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);

	//uv buffer
	cont.glVBO[v_quad_uv] = createBuffer();
	uploadBuffer(GL_ARRAY_BUFFER, &cont.uvs[0], cont.n_uvs * sizeof(float), cont.glVBO[v_quad_uv], false);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);

	//indices
	cont.glVBO[v_quad_ind] = createBuffer();
	uploadBuffer(GL_ELEMENT_ARRAY_BUFFER, &cont.inds[0], cont.n_inds * sizeof(unsigned int), cont.glVBO[v_quad_ind], false);

	glBindVertexArray(0);
}

//initialize pipeline
bool initFrameBuffer(Context& rend, Content& cont) {
	//The pathtrace frame buffer
	cont.glFBO[f_frame] = createFrameBuffer(int(rend.viewp[2]), int(rend.viewp[3]), cont.glRenTex[r_frame], cont.glRenTex[r_component], false);

	//best check ever
	GLuint status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT ||
		status == GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT ||
		status == GL_FRAMEBUFFER_UNSUPPORTED)
	{
		printf("Failed to generate context %i\n", status);
		return false;
	}

	//unbind frame buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	return true;
}

//do the path tracing
void renderContent(Context& rend, Content& cont) {
	//render to screen
	glBindFramebuffer(GL_FRAMEBUFFER, cont.glFBO[f_frame]);
	//clear the screen
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Use our shader?
	if (cont.shade == nullptr) {
		printf("no content shader\n");
		return;
	}
	Shader& shade = *cont.shade;
	glUseProgram(shade.getID());

	//use vao
	glBindVertexArray(cont.glVAO[a_quad]);

	//some uniforms...
	glUniform2f(shade.uniform("iResolution"), rend.viewp[2], rend.viewp[3]);
	glUniform1i(shade.uniform("iFrame"), rend.fnum);
	glUniform1f(shade.uniform("iGlobalTime"), rend.globalTime);
	glUniform3f(shade.uniform("loc"), loc.x, loc.y, loc.z);
	glUniform3f(shade.uniform("vel"), vel.x, vel.y, vel.z);
	glUniform3f(shade.uniform("iMouse"), float(mouseX), float(mouseY), float(mouseP));
	glUniform3f(shade.uniform("orient"), orient.x, orient.y, orient.z);

	//bind texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_frame]);
	glUniform1i(shade.uniform("iChannel0"), 0);
	//glActiveTexture(GL_TEXTURE1);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_depth]);
	//glUniform1i(shade.uniform("iChannel1"), 1);
	//glActiveTexture(GL_TEXTURE2);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_vert]);
	//glUniform1i(shade.uniform("iChannel2"), 2);
	//glActiveTexture(GL_TEXTURE3);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_ray]);
	//glUniform1i(shade.uniform("iChannel3"), 3);

	// Draw the indexed triangles!
	glDrawElements(GL_TRIANGLES, sizeof(quad_inds), GL_UNSIGNED_INT, 0);

	glBindVertexArray(0);
	glDisable(cont.glRenTex[r_frame]);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

//program entry point
int main(int argc, char* argv[]) {
	bool running = true;
	cv::Mat frame;

	printf("Renderer\n");
	Context renderer = Context(640, 480, 60.f, "Kyle's Path Tracer");
	renderer.ClearColor(vec4(0.f, 0.f, 0.f, 1.f));
	if (!renderer.initPassthrough())
		return 1;
	
	printf("Shaders\n");
	Shader pathtracer("passthrough.vert", "pathtrace.frag");
	std::vector<std::string> ptUniforms = { "iChannel0", "iChannel1", "iChannel2", "iChannel3", "iResolution", "iFrame", "iGlobalTime", "loc", "vel", "iMouse", "orient" };
	pathtracer.addUniforms(ptUniforms);

	printf("Content\n");
	Content cont;
	initGeometry(cont);
	cont.shade = &pathtracer;
	
	if (!initFrameBuffer(renderer, cont))
		return 2;

	printf("\nDone.\n\n");
	//ShowWindow(GetConsoleWindow(), SW_HIDE);
	
	do {
		renderer.fnum ++;
		renderer.globalTime += 0.017f;

		vec3 tmp = rotateY(vel, orient.xy);
		loc += tmp;

		//path trace
		renderContent(renderer, cont);
		//put restults on screen
		renderer.display(cont.glRenTex[r_frame]);

		//block until done hopefully
		glFinish();
		//get input
		glfwPollEvents();
		
		if (glfwGetMouseButton(renderer.window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
			//last vars
			double lastMX = mouseX;
			double lastMY = mouseY;
			//continuous mode
			glfwGetCursorPos(renderer.window, &mouseX, &mouseY);

			//if already pressed
			if (mouseP > 0.) {
				//update orient
				orient.y += float(mouseX - lastMX) * LOOK;
				orient.x += float(lastMY - mouseY) * LOOK;
				//overflow
				if (orient.x < -pi_5) orient.x = -pi_5;
				if (orient.x > pi_5) orient.x = pi_5;
				if (orient.y < -pi) orient.y += pi2;
				if (orient.y > pi) orient.y -= pi2;
			}
			//mouse pressed
			mouseP = 1.0;
			
		} else {
			mouseP = 0.0;
		}

		if (glfwGetKey(renderer.window, GLFW_KEY_UP) == GLFW_PRESS || glfwGetKey(renderer.window, GLFW_KEY_W) == GLFW_PRESS) {
			vel.z += ACCEL;
		} else if (glfwGetKey(renderer.window, GLFW_KEY_DOWN) == GLFW_PRESS || glfwGetKey(renderer.window, GLFW_KEY_S) == GLFW_PRESS) {
			vel.z -= ACCEL;
		} else if (abs(vel.z) >= ACCEL - sml) {
			vel.z -= sign(vel.z) * ACCEL;
		}
		if (glfwGetKey(renderer.window, GLFW_KEY_LEFT) == GLFW_PRESS || glfwGetKey(renderer.window, GLFW_KEY_A) == GLFW_PRESS) {
			vel.x -= ACCEL;
		} else if (glfwGetKey(renderer.window, GLFW_KEY_RIGHT) == GLFW_PRESS || glfwGetKey(renderer.window, GLFW_KEY_D) == GLFW_PRESS) {
			vel.x += ACCEL;
		} else if (abs(vel.x) >= ACCEL - sml) {
			vel.x -= sign(vel.x) * ACCEL;
		}
		if (glfwGetKey(renderer.window, GLFW_KEY_SPACE) == GLFW_PRESS) {
			vel.y += ACCEL;
		} else if (glfwGetKey(renderer.window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
			vel.y -= ACCEL;
		} else if (abs(vel.y) >= ACCEL - sml) {
			vel.y -= sign(vel.y) * ACCEL;
		}
		if (glfwGetKey(renderer.window, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwWindowShouldClose(renderer.window) != 0)
			running = false;

	} while (running);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}
