#version 330 core
precision highp int;
precision highp float;

////////////////////////////// Geometry Buffer ////////////////////////////// 

layout(location = 0) out vec4 fragColor;

smooth in vec2 uv;

uniform sampler2D iChannel0; // ----
uniform sampler2D iChannel1; // ----
uniform sampler2D iChannel2; // ----
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

mat3 cameraPoses(in float t) {
    float time = mod(t, 6.0);
    if (time < 2.)
        return mat3(vec3( 4.8, 0.5, -9.5), vec3(0.20, 0.85, 0.0), vec3(0.0));
    else if (time < 4.)
        return mat3(vec3( 4.8, 0.5, -4.8), vec3(0.15, 2.33, 0.0), vec3(0.0));
    else if (time < 6.)
        return mat3(vec3(-3.5, 2.5, -4.0), vec3(0.10, 1.80, 0.0), vec3(0.0));
}

vec4 renderG(in vec2 fragCoord, in vec2 iResolution, in int iFrame) {
    
    float asp = iResolution.x / iResolution.y;
    vec2 ndca = (2.0 * fragCoord.xy / iResolution.xy - 1.0) * vec2(asp, 1.0);

    // get current cam
	vec3 l = loc.xyz;
    vec2 o = orient.xy;

    // animated camera
    // if (iFrame > 1) {
    //    // just smoothstep between some nice angles
    //    float t = iTime * 0.5;
    //    mat3 cLast = cameraPoses(t);
    //    mat3 cNext = cameraPoses(t + 1.0);
    //    float ft = fract(t);
    //    ft = ft*ft*(3.0-2.0*ft);
    //    l = mix(cLast[0], cNext[0], ft);
    //    o = mix(cLast[1], cNext[1], ft);
    // }
    
    // save current camera info
    if (int(fragCoord.y) == int(iResolution.y) - 1) {
        int x = int(fragCoord.x);
        if (x == 0)
            return vec4(l, 0.0);
        else if (x == 1)
            return vec4(o, 0.0, 0.0);
    }
    
    // initial scene march
    vec3 d = rotateXY(normalize(vec3(ndca, FOV)), o);
    vec2 h = march(l, d, -1);
    // encode normal, object ID, and depth
    vec4 ret = vec4(1.0);
    encodeGbuffer(ret, norcurv(l + d * h.x, eps).xyz, int(h.y), h.x - eps);
    return ret;
}

void main() {
    fragColor = renderG(gl_FragCoord.xy, iResolution, iFrame);
}
