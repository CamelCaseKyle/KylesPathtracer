#version 330 core
precision highp int;
precision highp float;

////////////////////////////// Specular Buffer ////////////////////////////// 

layout(location = 0) out vec4 fragColor;

smooth in vec2 uv;

uniform sampler2D iChannel0; // g buffer
uniform sampler2D iChannel1; // last frame
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

vec4 renderSpecular(in vec2 fragCoord, in vec2 iResolution, in int iFrame, in sampler2D gBuffer, in sampler2D lastFrame) {
    bool biased = false;
#ifdef BIASED
    // force biased rendering
    biased = true;
#endif
    // defines variables for current ray, last ray, current hit, velocity, etc
    decodeAll(gBuffer, lastFrame);
    
    // save current camera info
    if (int(fragCoord.y) == int(iResolution.y) - 1) {
        int x = int(fragCoord.x);
        if (x == 0)
            return vec4(rl, 0.0);
        else if (x == 1)
            return vec4(ro, 0.0, 0.0);
    }

    // surface curvature to virtual image power
    vec4 nc = norcurv(hl, eps);
    float lightDistance = length(hl - light.xyz);
    vec3 sl = hl + rd * mix(0.0, lightDistance, eps / sqrt(max(eps, nc.w)));
    // reproject last frame onto this one
    vec4 fcol = reproject(ll, lo, sl, ho, asp, lastFrame, iResolution);
    fcol.a = floor(fcol.a);

    // temporal smoothing
    int lvv = min(TEMPORALSMOOTHING-1, int(float(TEMPORALSMOOTHING) * 2.0 * sqrt(length(vv))));
    if (fcol.a > float(TEMPORALSMOOTHING-lvv))
        fcol *= float(TEMPORALSMOOTHING-lvv) / fcol.a;

    // sample from surface
    mat3 surf = getSurface(ho, hl);
    // emissive
    fcol.rgb += surf[1];
    // not nothing
    if (ho != LIGHT) {
        int seed = genSeed(iFrame, ivec2(fragCoord), ivec2(iResolution));
        // do biased sampling
        if (biased) {
            // multiple importance sample lights and surfaces
            fcol.rgb += SMIS(rd, hl, hn, ho, seed);
        } else {
            // unbiased sampling (GT)
            fcol.rgb += UnbiasedPhong(rd, hl, hn, ho, seed);
        }
    }
#ifndef BIASED
    if (biased)
        fcol *= float(BIAS_WEIGHT);
#endif
    fcol.a += 1.0 + encodeBuffer(ho);
    return fcol;
}

void main() {
    fragColor = renderSpecular(gl_FragCoord.xy, iResolution, iFrame, iChannel0, iChannel1);
}
