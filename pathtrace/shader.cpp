#include "shader.h"

//some helper functions
std::string shaderFromFile(std::string file_path) {
	// Read the Shader code from the file
	std::string shaderCode;
	std::ifstream shaderStream(file_path, std::ios::in);
	if (shaderStream.is_open()){
		std::string Line = "";
		while (std::getline(shaderStream, Line))
			shaderCode += "\n" + Line;
		shaderStream.close();
	}
	else {
#ifdef DBG_SHADER
		printf("Impossible to open %s\n", file_path.c_str());
#endif
	}
	return shaderCode;
}

bool compileShader(GLuint shaderID, std::string& shaderCode) {
	GLint Result = GL_FALSE;
	int InfoLogLength;
	// Compile Geometry Shader
	char const *sourcePointer = shaderCode.c_str();
	glShaderSource(shaderID, 1, &sourcePointer, NULL);
	glCompileShader(shaderID);

	// Check Geometry Shader
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> shaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, NULL, &shaderErrorMessage[0]);
		return false;
	}
	return true;
}

//creates and loads a shader
Shader::Shader(std::string v, std::string f, bool feedback) {
	init = loadShader(v, f, feedback);
}
//adds a uniform to the pipeline
void Shader::addUniform(std::string n) {
	uniforms[n] = glGetUniformLocation(id, n.c_str());
}
//adds a few uniforms
void Shader::addUniforms(std::vector<std::string> n) {
	for (int i = 0; i < n.size(); i++) uniforms[n[i]] = glGetUniformLocation(id, n[i].c_str());
}
//loads a .vert .frag pair
bool Shader::loadShader(std::string vertex_file_path, std::string fragment_file_path, bool feedback) {
	bool win = true;

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if (VertexShaderStream.is_open()){
		std::string Line = "";
		while (std::getline(VertexShaderStream, Line))
			VertexShaderCode += "\n" + Line;
		VertexShaderStream.close();
	}
	else {
#ifdef DBG_SHADER
		printf("Impossible to open %s\n", vertex_file_path.c_str());
#endif
		//getchar();
		return false;
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if (FragmentShaderStream.is_open()){
		std::string Line = "";
		while (std::getline(FragmentShaderStream, Line))
			FragmentShaderCode += "\n" + Line;
		FragmentShaderStream.close();
	}


	GLint Result = GL_FALSE;
	int InfoLogLength;

	// Compile Vertex Shader
	char const *VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer, NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		//printf("%s\n", &VertexShaderErrorMessage[0]);
		win = false;
	}

	// Compile Fragment Shader
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer, NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		//printf("%s\n", &FragmentShaderErrorMessage[0]);
		win = false;
	}

	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);

	if (feedback) {
		const GLchar* feedbackVaryings[] = { "feedback" };
		glTransformFeedbackVaryings(ProgramID, 1, feedbackVaryings, GL_INTERLEAVED_ATTRIBS);
	}

	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);

#ifdef DBG_SHADER
		printf("%s / %s:\n%s\n", vertex_file_path.c_str(), fragment_file_path.c_str(), &ProgramErrorMessage[0]);
#endif
		win = false;
	}

	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

	id = ProgramID;

	init = true;
	return true;
}
//loads a .geo .vert .frag sequence
bool Shader::loadShader(std::string geometry_file_path, std::string vertex_file_path, std::string fragment_file_path, bool feedback) {
	// Create the shaders
	GLuint GeometryShaderID = glCreateShader(GL_GEOMETRY_SHADER);
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	//load the shaders
	std::string GeometryShaderCode = shaderFromFile(geometry_file_path);
	std::string VertexShaderCode = shaderFromFile(vertex_file_path);
	std::string FragmentShaderCode = shaderFromFile(fragment_file_path);
	//compile the shaders
	bool win = true;
	win &= compileShader(GeometryShaderID, GeometryShaderCode);
	win &= compileShader(VertexShaderID, VertexShaderCode);
	win &= compileShader(FragmentShaderID, FragmentShaderCode);

	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, GeometryShaderID);
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);

	if (feedback) {
		const GLchar* feedbackVaryings[] = { "feedback" };
		glTransformFeedbackVaryings(ProgramID, 1, feedbackVaryings, GL_INTERLEAVED_ATTRIBS);
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;

	glLinkProgram(ProgramID);
	//check program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);

#ifdef DBG_SHADER
		printf("Program error:\n%s\n", &ProgramErrorMessage[0]);
#endif
		win = false;
	}
	//clean up
	glDeleteShader(GeometryShaderID);
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

	id = ProgramID;

	init = true;
	return true;
}
