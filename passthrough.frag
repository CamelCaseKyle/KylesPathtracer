#version 330 core
precision highp int;
precision highp float;

////////////////////////////// Image Buffer ////////////////////////////// 

layout(location = 0) out vec4 fragColor;

smooth in vec2 uv;

uniform sampler2D iChannel0; // g buffer
uniform sampler2D iChannel1; // diffuse buffer
uniform sampler2D iChannel2; // specular buffer
uniform sampler2D iChannel3; // ----
uniform vec2 iResolution;
uniform int iFrame;
uniform float iTime;
uniform vec3 loc;
uniform vec3 vel;
uniform vec3 iMouse;
uniform vec3 orient;

#include "common.glsl"
// ignore this error with the linter

// exposure adjustment
const float brightness = 10.0;

vec4 renderImage(in vec2 fragCoord, in vec2 iResolution, in int iFrame, sampler2D gBuffer, sampler2D diffuse, sampler2D specular) {
    // defines variables for current ray, last ray, current hit, velocity, etc
    decodeAll(gBuffer, diffuse);

    // sample Diffuse buffer
    vec4 dBuffer = texelFetch(diffuse, ivec2(fragCoord), 0);
    // sample Specular buffer
    vec4 sBuffer = texelFetch(specular, ivec2(fragCoord), 0);
    
    // apply first surface properties
    mat3 surf = getSurface(ho, hl);
    dBuffer.rgb *= surf[0] * surf[2].x;
    sBuffer.rgb *= sqrt(surf[0]) * surf[2].y;
    
    // tonemap
    vec4 fcol = dBuffer / floor(dBuffer.a) + sBuffer / floor(sBuffer.a);
    fcol.rgb = linear_srgb(ACESFitted(fcol.rgb * brightness));
    return fcol;
}

void main() {
    fragColor = renderImage(gl_FragCoord.xy, iResolution, iFrame, iChannel0, iChannel1, iChannel2);
}
