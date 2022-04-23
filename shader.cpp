#include "shader.h"

using namespace std;

// trim from end of string
void rtrim(string &s, const char *t = " \t\f\v\n\r") {
	s.erase(s.find_last_not_of(t) + 1);
}

// trim from beginning of string (left)
void ltrim(string &s, const char *t = " \t\f\v\n\r") {
	s.erase(0, s.find_first_not_of(t));
}

// trim from both ends of string (left & right)
string& trim(string &s, const char *t = " \t\f\v\n\r") {
	rtrim(s, t);
	ltrim(s, t);
	return s;
}

// split a string around a character as a delimiter
vector<string> split(std::string &s, const char *t = " ") {
	vector<string> ret;

	size_t start = 0;
	size_t len = s.length();
	size_t lt = strlen(t);
	size_t end = s.find(t);
	while (end < len) {
		ret.push_back(s.substr(start, end - start));
		start = end + lt;
		end = s.find(t, start);
	}
	ret.push_back(s.substr(start, end));
	return ret;
}

// load shader and handle IO errors
string shaderFromFile(string file_path, bool recurse = true) {
	// Read the Shader code from the file
	string shaderCode = "";
	ifstream shaderStream(file_path, ios::in);
	if (shaderStream.is_open()){
		string Line = "", tLine = "";
		while (getline(shaderStream, Line)) {
			tLine = trim(Line);
			if (tLine.substr(0, 9) == "#include ") {
				vector<string> incLine = split(tLine);
				// copy contents into shader file (max 1 recur)
				if (recurse && incLine.size() > 1)
					shaderCode += "\n" + shaderFromFile(trim(incLine[1],"\""), false);
			} else {
				shaderCode += "\n" + Line;
			}
		}
		shaderStream.close();
	} else {
#ifdef DBG_SHADER
		printf("Cannot open %s\n", file_path.c_str());
#endif
	}
	return shaderCode;
}

// compile shader and output errors to std out
bool compileShader(string &shaderFile, GLuint shaderID, string &shaderCode) {
	GLint Result = GL_FALSE;
	int InfoLogLength;
	// Compile Shader
	char const *sourcePointer = shaderCode.c_str();
	glShaderSource(shaderID, 1, &sourcePointer, NULL);
	glCompileShader(shaderID);
	// Check Shader
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		vector<char> shaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, NULL, &shaderErrorMessage[0]);
		printf("Shader Compile Error: %s\n%s\n\n", shaderFile.c_str(), &shaderErrorMessage[0]);
		return false;
	}
	return true;
}

//creates and loads a shader
Shader::Shader(string v, string f, bool feedback) {
	init = loadShader(v, f, feedback);
}

//adds a uniform to the pipeline
void Shader::addUniform(string n) {
	uniforms[n] = glGetUniformLocation(id, n.c_str());
}

//adds a few uniforms
void Shader::addUniforms(vector<string> n) {
	for (int i = 0; i < n.size(); i++) uniforms[n[i]] = glGetUniformLocation(id, n[i].c_str());
}

//loads a .vert .frag pair
bool Shader::loadShader(string vertex_file_path, string fragment_file_path, bool feedback) {
	bool win = true;
	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	// Read the Shader code from the file
	string VertexShaderCode = shaderFromFile(vertex_file_path);
	string FragmentShaderCode = shaderFromFile(fragment_file_path);
	// Compile
	win &= compileShader(vertex_file_path, VertexShaderID, VertexShaderCode);
	win &= compileShader(fragment_file_path, FragmentShaderID, FragmentShaderCode);

	if (!win) return false;

	// create the program
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	if (feedback) {
		const GLchar* feedbackVaryings[] = { "feedback" };
		glTransformFeedbackVaryings(ProgramID, 1, feedbackVaryings, GL_INTERLEAVED_ATTRIBS);
	}
	glLinkProgram(ProgramID);

	// Check the program
	int InfoLogLength;
	GLint Result = GL_FALSE;
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
#ifdef DBG_SHADER
		printf("%s / %s:\n%s\n", vertex_file_path.c_str(), fragment_file_path.c_str(), &ProgramErrorMessage[0]);
#endif
		win = false;
	}
	// clean up
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);
	id = ProgramID;
	init = win;
	return win;
}

//loads a .geo .vert .frag sequence
bool Shader::loadShader(string geometry_file_path, string vertex_file_path, string fragment_file_path, bool feedback) {
	bool win = true;
	// Create the shaders
	GLuint GeometryShaderID = glCreateShader(GL_GEOMETRY_SHADER);
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	// load the shaders
	string GeometryShaderCode = shaderFromFile(geometry_file_path);
	string VertexShaderCode = shaderFromFile(vertex_file_path);
	string FragmentShaderCode = shaderFromFile(fragment_file_path);
	// compile the shaders
	win &= compileShader(geometry_file_path, GeometryShaderID, GeometryShaderCode);
	win &= compileShader(vertex_file_path, VertexShaderID, VertexShaderCode);
	win &= compileShader(fragment_file_path, FragmentShaderID, FragmentShaderCode);

	if (!win) return false;

	// create the program
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, GeometryShaderID);
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	if (feedback) {
		const GLchar* feedbackVaryings[] = { "feedback" };
		glTransformFeedbackVaryings(ProgramID, 1, feedbackVaryings, GL_INTERLEAVED_ATTRIBS);
	}
	glLinkProgram(ProgramID);

	// check program
	GLint Result = GL_FALSE;
	int InfoLogLength;
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
#ifdef DBG_SHADER
		printf("Program error:\n%s\n", &ProgramErrorMessage[0]);
#endif
		win = false;
	}
	// clean up
	glDeleteShader(GeometryShaderID);
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);
	id = ProgramID;
	init = win;
	return win;
}
