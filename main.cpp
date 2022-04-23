#define WIN32_LEAN_AND_MEAN

#include "render.h"

#include <Windows.h>
#include <chrono>

using namespace std;
using namespace chrono;
using namespace glm;

// index helpers
const enum fboIDs	  { f_gBuff, f_diffuse, f_specular };
const enum rentexIDs  { r_gBuff, r_diffuse, r_specular };
const enum shaderIDs  { s_gBuff, s_diffuse, s_specular };
const enum vaoIDs	  { a_quad };
const enum vboIDs	  { v_quad, v_quad_uv, v_quad_ind };

const float quad_verts[8] = {
	-1.0f, -1.0f,
	-1.0f, 1.0f,
	1.0f, -1.0f,
	1.0f, 1.0f,
};
const float quad_uvs[8] = {
	0.0f, 0.0f,
	0.0f, 1.0f,
	1.0f, 0.0f,
	1.0f, 1.0f,
};
const unsigned int quad_inds[4] = {
	0, 1, 2, 3
};

static const float
	accelSpeed = 0.01f, rotSpeed = 0.002f, maxSpeed = 0.5f,
	zfar = 1000.0f, eps = 0.00001f, ieps = 0.99999f, sml = 0.001f, isml = 0.999f,
	sc45 = 0.7071067f,
	HPI = 1.5707963f, PI2 = 6.2831853f, PI = 3.1415926f, iPI = 0.3183098f;

vec3 loc = vec3(-2.0f, 2.5f, -5.0f),
	 vel = vec3(0.0f, 0.0f, 0.0f),
	 orient = vec3(0.1f, 1.8f, 0.0f);
double mouseX, mouseY, mouseP;

inline float sign(float val) { return float((0.0f < val) - (val < 0.0f)); }

vec3 rotateY(const vec3 &p, const vec2 &angle) {
	vec2 c = cos(angle), s = sin(angle); 
	vec3 o = p;
	//o.yz = o.yz * mat2(c.x, s.x, -s.x, c.x);
	o.xz = o.xz * mat2(c.y, s.y, -s.y, c.y);
	return o;
}

// init content
void initGeometry(Content &cont) {
	cont.model = mat4(1.0);
	cont.view = mat4(1.0);
	cont.view_ray = mat4(1.0);
	cont.normalMatrix = mat4(1.0);
	cont.verts = quad_verts;
	cont.uvs = quad_uvs;
	cont.inds = quad_inds;
	cont.n_verts = 4 * 2;
	cont.n_uvs = 4 * 2;
	cont.n_inds = 4;

	cont.glVAO[a_quad] = createVertexArray(false);

	// vert buffer
	cont.glVBO[v_quad] = createBuffer();
	uploadBuffer(GL_ARRAY_BUFFER, &cont.verts[0], cont.n_verts * sizeof(float), cont.glVBO[v_quad], false);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);

	// uv buffer
	cont.glVBO[v_quad_uv] = createBuffer();
	uploadBuffer(GL_ARRAY_BUFFER, &cont.uvs[0], cont.n_uvs * sizeof(float), cont.glVBO[v_quad_uv], false);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);

	// indices
	cont.glVBO[v_quad_ind] = createBuffer();
	uploadBuffer(GL_ELEMENT_ARRAY_BUFFER, &cont.inds[0], cont.n_inds * sizeof(unsigned int), cont.glVBO[v_quad_ind], false);

	glBindVertexArray(0);
}

// initialize pipeline
bool initFrameBuffer(Context &rend, Content &cont) {
	// The pathtrace frame buffer
	cont.glFBO[f_gBuff] = createFrameBuffer(int(rend.viewp[2]), int(rend.viewp[3]), cont.glRenTex[r_gBuff], false); //  had depth attachment
	// The pathtrace frame buffer
	cont.glFBO[f_diffuse] = createFrameBuffer(int(rend.viewp[2]), int(rend.viewp[3]), cont.glRenTex[r_diffuse], false); //  had depth attachment
	// The pathtrace frame buffer
	cont.glFBO[f_specular] = createFrameBuffer(int(rend.viewp[2]), int(rend.viewp[3]), cont.glRenTex[r_specular], false); //  had depth attachment

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

// do the path tracing
void renderGbuffer(Context &rend, Content &cont, Shader &shade) {
	// render to G buffer
	glBindFramebuffer(GL_FRAMEBUFFER, cont.glFBO[f_gBuff]);
	// Use our shader
	glUseProgram(shade.getID());
	// use vao
	glBindVertexArray(cont.glVAO[a_quad]);
	// some uniforms...
	glUniform2f(shade.uniform("iResolution"), rend.viewp[2], rend.viewp[3]);
	glUniform1i(shade.uniform("iFrame"), rend.fnum);
	glUniform1f(shade.uniform("iTime"), rend.globalTime);
	glUniform3f(shade.uniform("loc"), loc.x, loc.y, loc.z);
	glUniform3f(shade.uniform("vel"), vel.x, vel.y, vel.z);
	glUniform3f(shade.uniform("iMouse"), float(mouseX), float(mouseY), float(mouseP));
	glUniform3f(shade.uniform("orient"), orient.x, orient.y, orient.z);
	// bind textures
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel0"), 0);
	//glActiveTexture(GL_TEXTURE1);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel1"), 1);
	//glActiveTexture(GL_TEXTURE2);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel2"), 2);
	//glActiveTexture(GL_TEXTURE3);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel3"), 3);

	// Draw the indexed triangles!
	glDrawElements(GL_TRIANGLE_STRIP, cont.n_inds, GL_UNSIGNED_INT, 0);

	//glDisable(GL_TEXTURE0);
	//glDisable(GL_TEXTURE1);
	//glDisable(GL_TEXTURE2);
	//glDisable(GL_TEXTURE3);
	//glDisable(cont.glFBO[f_gBuff]);
	//glDisable(cont.glRenTex[r_gBuff]);
}

// do the diffuse buffer
void renderDbuffer(Context &rend, Content &cont, Shader &shade) {
	// render to diffuse buffer
	glBindFramebuffer(GL_FRAMEBUFFER, cont.glFBO[f_diffuse]);
	// Use our shader
	glUseProgram(shade.getID());
	// use vao
	glBindVertexArray(cont.glVAO[a_quad]);
	// some uniforms...
	glUniform2f(shade.uniform("iResolution"), rend.viewp[2], rend.viewp[3]);
	glUniform1i(shade.uniform("iFrame"), rend.fnum);
	glUniform1f(shade.uniform("iTime"), rend.globalTime);
	glUniform3f(shade.uniform("loc"), loc.x, loc.y, loc.z);
	glUniform3f(shade.uniform("vel"), vel.x, vel.y, vel.z);
	glUniform3f(shade.uniform("iMouse"), float(mouseX), float(mouseY), float(mouseP));
	glUniform3f(shade.uniform("orient"), orient.x, orient.y, orient.z);
	// bind textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_gBuff]);
	glUniform1i(shade.uniform("iChannel0"), 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_diffuse]);
	glUniform1i(shade.uniform("iChannel1"), 1);
	//glActiveTexture(GL_TEXTURE2);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel2"), 2);
	//glActiveTexture(GL_TEXTURE3);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel3"), 3);

	// Draw the indexed triangles!
	glDrawElements(GL_TRIANGLE_STRIP, cont.n_inds, GL_UNSIGNED_INT, 0);

	//glDisable(GL_TEXTURE0);
	//glDisable(GL_TEXTURE1);
	//glDisable(GL_TEXTURE2);
	//glDisable(GL_TEXTURE3);
	//glDisable(cont.glFBO[f_gBuff]);
	//glDisable(cont.glRenTex[r_gBuff]);
	//glDisable(cont.glRenTex[r_diffuse]);
}

// do the specular buffer
void renderSbuffer(Context &rend, Content &cont, Shader &shade) {
	// render to specular buffer
	glBindFramebuffer(GL_FRAMEBUFFER, cont.glFBO[f_specular]);
	// Use our shader
	glUseProgram(shade.getID());
	// use vao
	glBindVertexArray(cont.glVAO[a_quad]);
	// some uniforms...
	glUniform2f(shade.uniform("iResolution"), rend.viewp[2], rend.viewp[3]);
	glUniform1i(shade.uniform("iFrame"), rend.fnum);
	glUniform1f(shade.uniform("iTime"), rend.globalTime);
	glUniform3f(shade.uniform("loc"), loc.x, loc.y, loc.z);
	glUniform3f(shade.uniform("vel"), vel.x, vel.y, vel.z);
	glUniform3f(shade.uniform("iMouse"), float(mouseX), float(mouseY), float(mouseP));
	glUniform3f(shade.uniform("orient"), orient.x, orient.y, orient.z);
	// bind textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_gBuff]);
	glUniform1i(shade.uniform("iChannel0"), 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, cont.glRenTex[r_specular]);
	glUniform1i(shade.uniform("iChannel1"), 1);
	//glActiveTexture(GL_TEXTURE2);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel2"), 2);
	//glActiveTexture(GL_TEXTURE3);
	//glBindTexture(GL_TEXTURE_2D, cont.glRenTex[]);
	//glUniform1i(shade.uniform("iChannel3"), 3);

	// Draw the indexed triangles!
	glDrawElements(GL_TRIANGLE_STRIP, cont.n_inds, GL_UNSIGNED_INT, 0);

	//glDisable(GL_TEXTURE0);
	//glDisable(GL_TEXTURE1);
	//glDisable(GL_TEXTURE2);
	//glDisable(GL_TEXTURE3);
	//glDisable(cont.glFBO[f_gBuff]);
	//glDisable(cont.glRenTex[r_gBuff]);
	//glDisable(cont.glRenTex[r_specular]);
}

void handleInput(Context &rend) {
	// handle mouse
	if (glfwGetMouseButton(rend.window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
		// last vars
		double lastMX = mouseX;
		double lastMY = mouseY;
		// continuous mode
		glfwGetCursorPos(rend.window, &mouseX, &mouseY);
		// if already pressed
		if (mouseP > 0.) {
			// update orient
			orient.y += float(mouseX - lastMX) * rotSpeed;
			orient.x += float(lastMY - mouseY) * rotSpeed;
			// overflow
			if (orient.x < -HPI) orient.x = -HPI;
			if (orient.x > HPI) orient.x = HPI;
			if (orient.y < -PI)  orient.y += PI2;
			if (orient.y > PI)  orient.y -= PI2;
		}
		// mouse pressed
		mouseP = 1.0;
	} else {
		mouseP = 0.0;
	}
	// handle keyboard
	vec3 newVel = vec3(0.0f);
	if (glfwGetKey(rend.window, GLFW_KEY_UP) == GLFW_PRESS           || glfwGetKey(rend.window, GLFW_KEY_W) == GLFW_PRESS) {
		newVel.z += 1.0f;
	} else if (glfwGetKey(rend.window, GLFW_KEY_DOWN) == GLFW_PRESS  || glfwGetKey(rend.window, GLFW_KEY_S) == GLFW_PRESS) {
		newVel.z -= 1.0f;
	}
	if (glfwGetKey(rend.window, GLFW_KEY_LEFT) == GLFW_PRESS         || glfwGetKey(rend.window, GLFW_KEY_A) == GLFW_PRESS) {
		newVel.x -= 1.0f;
	} else if (glfwGetKey(rend.window, GLFW_KEY_RIGHT) == GLFW_PRESS || glfwGetKey(rend.window, GLFW_KEY_D) == GLFW_PRESS) {
		newVel.x += 1.0f;
	}
	if (glfwGetKey(rend.window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		newVel.y += 1.0f;
	} else if (glfwGetKey(rend.window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
		newVel.y -= 1.0f;
	}
	if (length(newVel) > 1.0f)
		newVel = normalize(newVel);
	// friction
	vel *= 0.9;
	// add velocity in our coord system
	vec3 tmp = rotateY(newVel*accelSpeed, orient.xy);
	vel += tmp;
	float speed = length(vel);
	// slow to complete stop
	if (speed < accelSpeed)
		vel *= 0.0f;
	// limit speed
	if (speed > maxSpeed)
		vel *= maxSpeed / speed;
	loc += vel;
}

// program entry point
int main(int argc, char* argv[]) {
	bool running = true;
	
	printf("Renderer\n");
	Context renderer = Context(1280, 720, 60.0f, "Kyle's Path Tracer");
	renderer.ClearColor(vec4(0.0f, 0.0f, 0.0f, 1.0f));
	if (!renderer.initPassthrough())
		return 1;
	
	printf("Shaders\n");
	vector<string> shadertoyUniforms = { "iChannel0", "iChannel1", "iChannel2", "iChannel3", "iResolution", "iFrame", "iGlobalTime", "loc", "vel", "iMouse", "orient" };
	Shader diffuseShader("passthrough.vert", "diffuse.frag");
	diffuseShader.addUniforms(shadertoyUniforms);
	Shader specularShader("passthrough.vert", "specular.frag");
	specularShader.addUniforms(shadertoyUniforms);
	Shader geometryShader("passthrough.vert", "geometry.frag");
	geometryShader.addUniforms(shadertoyUniforms);

	printf("Content\n");
	Content cont;
	initGeometry(cont);
	
	if (!initFrameBuffer(renderer, cont))
		return 2;

	printf("\nDone.\n\n");
	//ShowWindow(GetConsoleWindow(), SW_HIDE);
	
	high_resolution_clock::time_point start = high_resolution_clock::now();

	do {
		// get input
		glfwPollEvents();

		// exit game
		if (glfwGetKey(renderer.window, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwWindowShouldClose(renderer.window) != 0)
			running = false;

		// game controls
		handleInput(renderer);

		// render next frame
		renderer.fnum++;
		renderer.globalTime = float((high_resolution_clock::now() - start).count());

		// g buffer
		renderGbuffer(renderer, cont, geometryShader);
		// d buffer
		renderDbuffer(renderer, cont, diffuseShader);
		// s buffer
		renderSbuffer(renderer, cont, specularShader);
		// put restults on screen
		renderer.display(cont.glRenTex[r_gBuff], cont.glRenTex[r_diffuse], cont.glRenTex[r_specular]);
		// swap after everything is finished (hopefully driver wont just ignore)
		glfwSwapBuffers(renderer.window);

		// 60fps, calls above are non blocking
		Sleep(16);

	} while (running);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}
