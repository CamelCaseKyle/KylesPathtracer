#pragma once

#define DBG_SHADER

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

//#include <glm.hpp>
//#include <gtc/matrix_transform.hpp>
//#include <gtx/quaternion.hpp>

//#include <opencv2/core.hpp>
//#include <opencv2/highgui.hpp>

class Shader {
	GLuint id = NULL;
	//use dictionary or something
	std::unordered_map<std::string, GLuint> uniforms;

public:
	bool init = false;

	Shader() {}
	//creates and loads a shader
	Shader(std::string v, std::string f, bool feedback = false);

	//loads a .vert .frag pair
	bool loadShader(std::string vertex_file_path, std::string fragment_file_path, bool feedback = false);
	//loads a .geo .vert .frag series
	bool loadShader(std::string geometry_file_path, std::string vertex_file_path, std::string fragment_file_path, bool feedback = false);

	//getters
	GLuint getID() { return id; }
	GLuint uniform(std::string name) { return uniforms[name]; }

	//adds a uniform to the pipeline
	void addUniform(std::string n);
	//add a few uniforms
	void addUniforms(std::vector<std::string> n);
};
