#version 330 core

precision highp int;
precision highp float;

layout(location = 0) out vec4 fragColor;

smooth in vec2 uv;

uniform sampler2D iChannel0;
uniform sampler2D iChannel1; 
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform vec2 iResolution;
uniform int iFrame;
uniform float iGlobalTime;
uniform vec3 loc;
uniform vec3 vel;
uniform vec3 iMouse;
uniform vec3 orient;

// Thanks Paniq
vec3 linear_srgb(vec3 x) {
    return mix(1.055*pow(x, vec3(1./2.4)) - 0.055, 12.92*x, step(x, vec3(0.0031308)));
}

vec3 srgb_linear(vec3 x) {
    return mix(pow((x + 0.055)/1.055,vec3(2.4)), x / 12.92, step(x, vec3(0.04045)));
}

// Paniq's ACES fitted from https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
vec3 ACESFitted(vec3 color) {
	// ODT_SAT => XYZ => D60_2_D65 => sRGB
    color = color * mat3(
        0.59719, 0.35458, 0.04823,
        0.07600, 0.90834, 0.01566,
        0.02840, 0.13383, 0.83777
    );
    // Apply RRT and ODT
    vec3 a = color * (color + 0.0245786) - 0.000090537;
    vec3 b = color * (0.983729 * color + 0.4329510) + 0.238081;
    color = a / b;
	// Back to color space
    color = color * mat3(
         1.60475, -0.53108, -0.07367,
        -0.10208,  1.10813, -0.00605,
        -0.00327, -0.07276,  1.07602
    );
    // Clamp to [0, 1]
    return clamp(color, 0.0, 1.0);
}

void main() {
    fragColor = texture2D(iChannel0, uv);
    fragColor.rgb = linear_srgb(ACESFitted(fragColor.rgb / floor(fragColor.a)));
}
