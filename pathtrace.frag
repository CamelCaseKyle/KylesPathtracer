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


//////////////////////////////// Rendering Configuration ////////////////////////////////

// always use biased sampling (fallback to unbiased ground truth)
#define BIASED
// use low discrepency or psuedo-random sampling
//#define STRUCTURED

// ray bounces
#define BOUNCES 1
// raymarching steps
#define STEPS 128
// number of frames to reproject and smooth over time
#define TEMPORALSMOOTHING 30

// biased light samples
#define SMP_DIRECT 1
// unbiased light samples
#define SMP_UNBIAS 8
// specular samples
#define SMP_IDS 1
// diffuse samples
#define SMP_IDD 1
// bias sample weight to reduce initial unbiased variance
#define BIAS_WEIGHT 1

//////////////////////////////// Math Toolkit ////////////////////////////////

const float	eps = 0.001,     ieps = 0.999,   zfar = 50.0,       FOV = 1.5,
            HPI = 1.5707963, PI = 3.1415926, TWOPI = 6.2831853, SQRT2 = 1.4142136, SC45 = 0.7071068;
const vec2	  VEL = vec2(0.5, 0.5),   POS = vec2(1.5, 0.5),   ROT = vec2(2.5, 0.5),   MOU = vec2(3.5, 0.5);
const vec2	L_VEL = vec2(4.5, 0.5), L_POS = vec2(5.5, 0.5), L_ROT = vec2(6.5, 0.5), L_MOU = vec2(7.5, 0.5);

// generate a unique value for each pixel/sample/frame
ivec3 seed(in int iFrame, in int smp, in ivec2 fragCoord) {
    ivec2 s = fragCoord + iFrame + smp;
    return ivec3(s, s.x ^ s.y);
}

float weyl1(in int v) {
    return fract(float(v*10368889) / exp2(24.0));
}
// thanks jaybird https://www.shadertoy.com/view/4dtBWH
vec2 weyl2(in int v) {
    return fract(vec2(v*ivec2(12664745, 9560333)) / exp2(24.0));
}
vec3 weyl3(in int v) {
    return fract(vec3(v*ivec3(13743434, 11258243, 9222443)) / exp2(24.0));
}

vec3 pcg3(in ivec3 v) {
    uvec3 w = uvec3(v) * 1664525u + 1013904223u;
    w.x += w.y*w.z;
    w.y += w.z*w.x;
    w.z += w.x*w.y;
    w ^= w >> 16u;
    w.x += w.y*w.z;
    w.y += w.z*w.x;
    w.z += w.x*w.y;
    return vec3(w) / exp2(32.0);
}

// fudged a bit to cut it off as close to the (0,1) interval as possible
vec3 logit3(in vec3 v) {
    vec3 t = 0.988 * (v + 0.006);
    return log(t / (1.0 - t)) * 0.221 + 0.5;
}

// thanks hornet https://www.shadertoy.com/view/4ssXRX
vec3 boxmuller(in vec3 u, in vec3 v) {
	return 0.23 * sqrt(-log(u + 0.00001))*cos(TWOPI * v) + 0.5;
}

void basis(in vec3 n, out vec3 f, out vec3 r) {
    float s = (n.z >= 0.0)? 1.0: -1.0;
    float a = 1.0 / (s + n.z);
    float b = -n.x*n.y*a;
    f = vec3(1.0 - n.x*n.x*a*s, b*s, -n.x*s);
    r = vec3(b, s - n.y*n.y*a, -n.y);
}

vec3 rotateXY(in vec3 p, in vec2 angle) {
	vec2 c = cos(angle), s = sin(angle);
    vec3 o = p;
	o.yz *= mat2(c.x, s.x, -s.x, c.x); 
    o.xz *= mat2(c.y, s.y, -s.y, c.y);
	return o;
}

// thanks iq http://iquilezles.org/www/articles/texture/texture.htm
vec4 textureGood(sampler2D sam, in vec2 x, in int bits) {
	ivec2 p = ivec2(floor(x));
    vec2 f = fract(x);
    f = f*f*(3.0-2.0*f);
    vec4 a = texelFetch(sam, (p+ivec2(0,0)) & bits, 0);
	vec4 b = texelFetch(sam, (p+ivec2(1,0)) & bits, 0);
	vec4 c = texelFetch(sam, (p+ivec2(0,1)) & bits, 0);
	vec4 d = texelFetch(sam, (p+ivec2(1,1)) & bits, 0);
	return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}

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

float inRange(float a, float x, float b) {
    return step(a, x) * step(x, b);
}

//input is normalized visible spectrum wavelength (0. = 400nm violet/uv, 1. = 700nm red/ir)
vec3 spectrum(in float x) {
    float t, l = x*300.+400.;
    float[] fr1 = float[] (400.,410.,545.,595.,650.,415.,475.,585.,400.,475.),
			fr2 = float[] (410.,475.,595.,650.,700.,475.,585.,639.,475.,560.),
        	dv1 = float[] (10., 65., 50., 55., 50., 60., 115.,54., 75., 85.);
    vec3[] c = vec3[] (vec3(0.,.33,-.2),vec3(.14,0.,-.13),vec3(0.,1.98,-1.),
        vec3(.98,.06,-.4),vec3(.65,-.84,.2),vec3(0.,0.,.8), vec3(.8,.76,-.8),
        vec3(.84,-.84,0.),vec3(0.,2.2,-1.5),vec3(.7,-1.,.3));
    vec3 r = vec3(0.);
    for (int i = 0; i < 5; i++) {
        t = (l-fr1[i]) / dv1[i];
        r.r += inRange(fr1[i],l,fr2[i])*(c[i].x + c[i].y*t + c[i].z*t*t);
    }
    for (int i = 5; i < 8; i++) {
        t = (l-fr1[i]) / dv1[i];
        r.g += inRange(fr1[i],l,fr2[i])*(c[i].x + c[i].y*t + c[i].z*t*t);
    }
    for (int i = 8; i < 10; i++) {
        t = (l-fr1[i]) / dv1[i];
        r.b += inRange(fr1[i],l,fr2[i])*(c[i].x + c[i].y*t + c[i].z*t*t);
    }
    return r*r;
}

// linear angle of sphere at distance D with radius R
float linearAngle(float d, float r) {
    return asin(clamp(r/d, eps, ieps));
}

// solid angle of sphere given distance squared and radius squared
float solidAngle(float d2, float r2) {
    return (1.0 - sqrt(1.0 - clamp(r2/d2, 0.0, 1.0))) * TWOPI;
}

float Schlick(in float r1, in float r2, in float vn) {
	float r0 = (r1 - r2) / (r1 + r2);
	return mix(r0*r0, 1.0, pow(1.0 - vn, 5.0));
}

float LambertianPDF(in vec3 hn, in vec3 nlv) {
    return 1.0 / max(eps, dot(nlv, hn));
}

// using normalized gaussian distribution
vec3 uniformDir(ivec3 seed) {
#ifdef STRUCTURED
    vec3 rnd = logit3(weyl3(seed.z));
#else
    vec3 rnd = boxmuller(pcg3(seed), pcg3(~seed.zxy));
#endif
    return normalize(rnd * 2.0 - 1.0);
}

// only one hemisphere
vec3 uniformHemiDir(vec3 hn, ivec3 seed) {
    vec3 rnd = uniformDir(seed);
    return rnd * sign(dot(hn, rnd));
}
float uniformHemiPDF() {
    return 1.0 / TWOPI;
}

// cosine distribution
vec3 cosHemiDir(vec3 hn, ivec3 seed) {
    vec3 rnd = uniformDir(seed);
    return normalize(hn + rnd * ieps);
}
float cosHemiPDF(in vec3 hn, in vec3 lv) {
    return max(eps, dot(lv, hn)) / PI;
}

// uniform sample cone
vec3 uniformConeDir(vec3 lv, float llv, float lr, in ivec3 seed) {
#ifdef STRUCTURED
    vec3 rnd = weyl3(seed.z);
#else
    vec3 rnd = pcg3(seed);
#endif
    vec3 nlv = lv / llv;
    float rad = sqrt(rnd.x);
    float tha = rnd.y * TWOPI;
    vec3 r, u;
    basis(nlv, r, u);
    r *= cos(tha);
    u *= sin(tha);
    float sa = linearAngle(llv, lr);
    return normalize(nlv + tan(sa) * rad * (r + u));
}
float uniformConePDF(in vec3 lv, in float lr) {
    return 1.0 / solidAngle(dot(lv, lv), lr*lr);
}

//////////////////////////////// Scene Modeling ////////////////////////////////

const int
    LIGHT = 1,
    FLOOR = 2,
    WALL1 = 3,
    WALL2 = 6,
    CEIL = 7,
    BOX = 8;

vec3 lightLoc = vec3(4., 5., -4.);
float lightRadius = 1.;
vec3 lightColor = vec3(100.);

mat3 getSurface(in int ho, in vec3 hl) {
    mat3 ret;
    if (ho == LIGHT) {
        ret[0] = vec3(0.0); // reflection color
        ret[1] = lightColor; // emission color
        ret[2] = vec3(0.5, 0.5, 0.0); // diffuse (x) specular (y)
    } else if (ho == BOX) {
        ret[0] = vec3(0.05 + 0.25 * float(int(floor(hl.x * 4.0) + floor(hl.y * 4.0) + floor(hl.z * 4.0)) & 1)); // reflection color
        ret[1] = vec3(0.0); // emission color
        ret[2] = vec3(0.5, 0.5, 0.0); // roughness (x) specular (y)
    } else if (ho <= 0) {

    } else {
        float refl = float(int(ho == FLOOR || ho == CEIL)) * (0.5 + float(int(floor(hl.x) + floor(hl.y) + floor(hl.z)) & 1)) * 0.25 + 0.75;
        float matCol = float(ho);
        ret[0] = vec3(0.05 + cos(matCol) * 0.025, 0.05 + sin(matCol) * 0.025, 0.05) * refl;
        ret[1] = vec3(0.0); // emission color
        ret[2] = vec3(refl, refl, 0.0); // roughness (x) specular (y)
    }
    return ret;
}

vec2 sdMin(vec2 a, vec2 b) {
    if (a.x < b.x)
        return a;
    return b;
}

float sdBox(in vec3 p, in mat3 o, in vec3 s) {
    vec3 d = abs(p * o) - s;
    return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
}

vec2 sdf(vec3 l) {
    vec2 d = vec2(zfar, 0.);
    
    // floor
    d = sdMin(d, vec2(l.y - 0.0, FLOOR));
    
    // ceil
    d = sdMin(d, vec2(-l.y + 10.0, CEIL));
    
    // walls
    d = sdMin(d, vec2(-l.x + 10.0, WALL1));
    d = sdMin(d, vec2(l.z + 10.0, WALL2));

    // rounded cube
    d = sdMin(d, vec2(sdBox(l - vec3(8.0, 1.0, -8.0), mat3(1.0), vec3(0.8)) - 0.1, BOX));
    d = sdMin(d, vec2(distance(l, vec3(8.0, 2.5, -8.0)) - 0.49, BOX));
    
    // light
    d = sdMin(d, vec2(length(l - lightLoc) - lightRadius, LIGHT));
    
    return d;
}

vec3 norm(in vec3 p, in float ep) {
    vec3 n = vec3(0.0);
    for (int i = 0; i < 4; ++i) {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1), ((i>>1)&1), (i&1))-1.0);
        n += e * sdf(p + e*ep).x;
    }
    return normalize(n);
}

vec3 march(in vec3 l, in vec3 rd) {
    float t = 0.0;
    vec2 sdSmp;
    for (int i = 0; i < STEPS; ++i) {
        sdSmp = sdf(l + rd * t);
        t += sdSmp.x * 0.9;
        if (sdSmp.x < eps)
            break;
        if (t > zfar)
            return vec3(zfar, 0.0, 0.0);
    }
    return vec3(min(t, zfar), 0.0, sdSmp.y);
}

vec3 lightMarch(in vec3 hl, in vec3 hn, in vec3 ld) {
    int ho = int(march(hl + hn * eps, ld).z);
    if (ho == LIGHT) {
        // hit the light
        return lightColor;
    } else {
        // hit nothing
        return vec3(0.0);
    }
}

//////////////////////////////// approx PDF only ////////////////////////////////

float lightPDF(out vec3 directDir, in vec3 hl, in vec3 hn, in vec3 ll, in float lr, in ivec3 seed) {
    vec3 lv = ll - hl;
    directDir = uniformConeDir(lv, length(lv), lr, seed);
    float lpdf = uniformConePDF(lv, lr);
    float gpdf = LambertianPDF(hn, directDir);
    return 1.0 / (lpdf * gpdf);
}

// samples a specular reflection off a plane
float specularPlanePDF(out vec3 specDir, in vec3 hl, in vec3 hn, in vec3 ll, in float lr, in vec3 pl, in vec3 pn, in ivec3 seed) {
    // get distance from each point to the plane
    float a = dot(hl - pl, pn),
          b = dot(ll - pl, pn);
    // similar triangles, just use the ratio of side lengths
    vec3 s = mix(hl - a*pn, ll - b*pn, a / (a + b));
    // surface to specular point s
    vec3 sv = s - hl;
    float sv2 = dot(sv, sv);
    float lsv = sqrt(sv2);
    // specular point to light
    vec3 ls = ll - s;
    float ls2 = dot(ls, ls);
    float lls = sqrt(ls2);
    // sample the image of the specular reflection on the surface TODO: implement roughness here
    vec3 ts = sv * lls;
    float ts2 = dot(ts, ts);
    float lts = sqrt(ts2);
    specDir = uniformConeDir(ts, lts, lr*lsv, seed);
    // light solid angle and schlick PDF
    float gpdf = LambertianPDF(hn, specDir);
    float lpdf = uniformConePDF(ts, lr*lsv);
    float spdf = 1.0 / Schlick(1.0003, 2.4, max(eps, dot(specDir, pn)));
    return 1.0 / (gpdf * lpdf * spdf);
}

float diffusePlanePDF(out vec3 diffDir, in vec3 hl, in vec3 hn, in vec3 ll, in float lr, in vec3 pl, in vec3 pn, in ivec3 seed) {
    // dot project light onto plane
    vec3 d = ll - pn * dot(ll - pl, pn);
    // surface to diffuse point d
    vec3 dv = d - hl;
    float dv2 = dot(dv, dv);
    float ldv = sqrt(dv2);
    // diffuse point to light
    vec3 ld = ll - d;
    float ld2 = dot(ld, ld);
    float lld = sqrt(ld2);
    // cosine sample disc on ground below
    diffDir = uniformConeDir(dv, ldv, min(ldv, lld) * 0.9, seed);
    float dpdf = PI;
    float g2pdf = LambertianPDF(pn, -diffDir);
    float gpdf = LambertianPDF(hn, diffDir);
    return 1.0 / (dpdf * gpdf * g2pdf);
}

//////////////////////////////// full PDF and sampling ////////////////////////////////

vec3 lightContribution(in vec3 directDir, in float density, in vec3 hl, in vec3 hn, in vec3 ll, in float lr, in ivec3 seed) {
    return lightMarch(hl, hn, directDir) * density;
}

vec3 specularPlaneContribution(in vec3 specDir, in float density, in vec3 hl, in vec3 hn, in vec3 pl, in vec3 pn, in int po, in vec3 ll, in float lr, in ivec3 seed) {
    // specular plane sample
    vec3 sres = march(hl + hn * eps, specDir);
    if (int(sres.z) != po)
        return vec3(0.0);
    vec3 _hl = hl + specDir * sres.x;
    vec3 _lv = lightLoc - _hl;
    float _lv2 = dot(_lv, _lv);
    // light sampling
    vec3 sampleDir = uniformConeDir(_lv, sqrt(_lv2), lightRadius, seed);
    vec3 li = lightMarch(_hl, pn, sampleDir);
    mat3 surf = getSurface(int(sres.z), _hl);
    // surface emission, color, roughness
    return surf[1] + surf[2].y * surf[0] * li * density;
}

vec3 diffusePlaneContribution(in vec3 diffDir, in float density, in vec3 hl, in vec3 hn, in vec3 pl, in vec3 pn, in int po, in vec3 ll, in float lr, in ivec3 seed) {
    // diffuse plane sample
    vec3 dres = march(hl + hn * eps, diffDir);
    if (int(dres.z) != po)
        return vec3(0.0);
    vec3 _hl = hl + diffDir * dres.x;
    vec3 _lv = lightLoc - _hl;
    float _lv2 = dot(_lv, _lv);
    // light sampling
    vec3 sampleDir = uniformConeDir(_lv, sqrt(_lv2), lightRadius, seed);
    vec3 li = lightMarch(_hl, pn, sampleDir);
    mat3 surf = getSurface(int(dres.z), _hl);
    // surface emission, color, roughness
    return surf[1] + surf[2].x * surf[0] * li * density;
}

//////////////////////////////// Sampling Strategy ////////////////////////////////

vec3 MIS(in vec3 hl, in vec3 hn, in int ho, in ivec3 s) {
    vec3 ret = vec3(0.0);
    
#if SMP_DIRECT
    // take direct light samples
    vec3 smpDirect = vec3(0.0);
    // get all PDFs
    vec4 dlpdf = vec4(0.0);
    if (ho != LIGHT)
        dlpdf.w = lightPDF(dlpdf.xyz, hl, hn, lightLoc, lightRadius, s);
    // compute CDF
    float lightCDF = dlpdf.w;
    // MIS weights
    float dlweight = lightCDF / max(eps, dlpdf.w);

    for (int i = 0; i < SMP_DIRECT; ++i) {
        ivec3 si = s + i;
  #ifdef STRUCTURED
        vec3 rnd = weyl3(si.z) * lightCDF;
  #else
        vec3 rnd = pcg3(si) * lightCDF;
  #endif
        if (rnd.z < lightCDF) {
            smpDirect += lightContribution(dlpdf.xyz, dlpdf.w, hl, hn, lightLoc, lightRadius, si);
        } else {
            smpDirect += vec3(1000.0, 0.0, 0.0);
        }
    }
    ret += smpDirect / float(SMP_DIRECT);
#endif

#if SMP_IDS
    // take indirect specular samples
    vec3 smpIDS = vec3(0.0);
    // get all PDFs
    vec4 fspdf = vec4(0.0);
    vec4 cspdf = vec4(0.0);
    vec4 w1spdf = vec4(0.0);
    vec4 w2spdf = vec4(0.0);
    if (ho != FLOOR)
        fspdf.w = specularPlanePDF(fspdf.xyz, hl, hn, lightLoc, lightRadius, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), s);
    if (ho != CEIL)
        cspdf.w = specularPlanePDF(cspdf.xyz, hl, hn, lightLoc, lightRadius, vec3(0.0, 10.0, 0.0), vec3(0.0, -1.0, 0.0), s);
    if (ho != WALL1)
        w1spdf.w = specularPlanePDF(w1spdf.xyz, hl, hn, lightLoc, lightRadius, vec3(10.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), s);
    if (ho != WALL2)
        w2spdf.w = specularPlanePDF(w2spdf.xyz, hl, hn, lightLoc, lightRadius, vec3(0.0, 0.0,-10.0), vec3(0.0, 0.0, 1.0), s);
    // compute CDF
    float floorSpecCDF = fspdf.w;
    float ceilSpecCDF = floorSpecCDF + cspdf.w;
    float w1SpecCDF = ceilSpecCDF + w1spdf.w;
    float w2SpecCDF = w1SpecCDF + w2spdf.w;
    // MIS weights outside for loop
    float fsweight = w2SpecCDF / max(eps, fspdf.w);
    float csweight = w2SpecCDF / max(eps, cspdf.w);
    float w1sweight = w2SpecCDF / max(eps, w1spdf.w);
    float w2sweight = w2SpecCDF / max(eps, w2spdf.w);

    for (int i = 0; i < SMP_IDS; ++i) {
        ivec3 si = s + i;
  #ifdef STRUCTURED
        vec3 rnd = weyl3(si.z) * w2SpecCDF;
  #else
        vec3 rnd = pcg3(si) * w2SpecCDF;
  #endif
        if (rnd.z <= floorSpecCDF) {
            smpIDS += specularPlaneContribution(fspdf.xyz, fspdf.w, hl, hn, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), FLOOR, lightLoc, lightRadius, si) * fsweight;
        } else if (rnd.z <= ceilSpecCDF) {
            smpIDS += specularPlaneContribution(cspdf.xyz, cspdf.w, hl, hn, vec3(0.0, 10.0, 0.0), vec3(0.0, -1.0, 0.0), CEIL, lightLoc, lightRadius, si) * csweight;
        } else if (rnd.z <= w1SpecCDF) {
            smpIDS += specularPlaneContribution(w1spdf.xyz, w1spdf.w, hl, hn, vec3(10.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), WALL1, lightLoc, lightRadius, si) * w1sweight;
        } else if (rnd.z <= w2SpecCDF) {
            smpIDS += specularPlaneContribution(w2spdf.xyz, w2spdf.w, hl, hn, vec3(0.0, 0.0, -10.0), vec3(0.0, 0.0, 1.0), WALL2, lightLoc, lightRadius, si) * w2sweight;
        } else {
            smpIDS += vec3(1000.0, 0.0, 0.0);
        }
    }
    ret += smpIDS / float(SMP_IDS);
#endif

#if SMP_IDD
    // take indirect diffuse samples
    vec3 smpIDD = vec3(0.0);
    // get all PDFs
    vec4 fdpdf = vec4(0.0);
    vec4 cdpdf = vec4(0.0);
    vec4 w1dpdf = vec4(0.0);
    vec4 w2dpdf = vec4(0.0);
    if (ho != FLOOR)
        fdpdf.w = diffusePlanePDF(fdpdf.xyz, hl, hn, lightLoc, lightRadius, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), s);
    if (ho != CEIL)
        cdpdf.w = diffusePlanePDF(cdpdf.xyz, hl, hn, lightLoc, lightRadius, vec3(0.0, 10.0, 0.0), vec3(0.0, -1.0, 0.0), s);
    if (ho != WALL1)
        w1dpdf.w = diffusePlanePDF(w1dpdf.xyz, hl, hn, lightLoc, lightRadius, vec3(10.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), s);
    if (ho != WALL2)
        w2dpdf.w = diffusePlanePDF(w2dpdf.xyz, hl, hn, lightLoc, lightRadius, vec3(0.0, 0.0, -10.0), vec3(0.0, 0.0, 1.0), s);
    // compute CDF
    float floorDiffCDF = fdpdf.w;
    float ceilDiffCDF = floorDiffCDF + cdpdf.w;
    float w1DiffCDF = ceilDiffCDF + w1dpdf.w;
    float w2DiffCDF = w1DiffCDF + w2dpdf.w;
    // MIS weights outside for loop
    float fdweight = w2DiffCDF / max(eps, fdpdf.w);
    float cdweight = w2DiffCDF / max(eps, cdpdf.w);
    float w1dweight = w2DiffCDF / max(eps, w1dpdf.w);
    float w2dweight = w2DiffCDF / max(eps, w2dpdf.w);
    
    for (int i = 0; i < SMP_IDD; ++i) {
        ivec3 si = s + i;
  #ifdef STRUCTURED
        vec3 rnd = weyl3(si.z) * w2DiffCDF;
  #else
        vec3 rnd = pcg3(si) * w2DiffCDF;
  #endif
        if (rnd.z <= floorDiffCDF) {
            smpIDD += diffusePlaneContribution(fdpdf.xyz, fdpdf.w, hl, hn, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), FLOOR, lightLoc, lightRadius, si) * fdweight;
        } else if (rnd.z <= ceilDiffCDF) {
            smpIDD += diffusePlaneContribution(cdpdf.xyz, cdpdf.w, hl, hn, vec3(0.0, 10.0, 0.0), vec3(0.0, -1.0, 0.0), CEIL, lightLoc, lightRadius, si) * cdweight;
        } else if (rnd.z <= w1DiffCDF) {
            smpIDD += diffusePlaneContribution(w1dpdf.xyz, w1dpdf.w, hl, hn, vec3(10.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), WALL1, lightLoc, lightRadius, si) * w1dweight;
        } else if (rnd.z <= w2DiffCDF) {
            smpIDD += diffusePlaneContribution(w2dpdf.xyz, w2dpdf.w, hl, hn, vec3(0.0, 0.0, -10.0), vec3(0.0, 0.0, 1.0), WALL2, lightLoc, lightRadius, si) * w2dweight;
        } else {
            smpIDD += vec3(1000.0, 0.0, 0.0);
        }
    }
    ret += smpIDD / float(SMP_IDD);
#endif
    
    return ret;
}

// ray reflect scatter
void brdf(inout vec3 ro, inout vec3 rd, in vec3 hl, in vec3 hn, in ivec3 seed) {
    ro = hl + hn * eps;
    rd = cosHemiDir(hn, seed);
}

//////////////////////////////// scene rendering ////////////////////////////////

float encode(in int o) {
    return float(o) * 0.101;
}

int decode(in float a) {
    return int(fract(a) * 10.10);
}

// last loc, last orient, hit loc, hit object, 1/aspect, channel, resolution
vec4 reproject(in vec3 ll, in vec3 lo, in vec3 hl, in int ho, in float asp, sampler2D iChannel, in vec2 iChanRes) {
    // last camera basis
    vec3 lf = rotateXY(vec3(0.0, 0.0, 1.0), lo.xy);
    vec3 r = normalize(cross(lf, vec3(0.0, 1.0, 0.0)));
    vec3 u = normalize(cross(lf, r));
    // dir to point
    vec3 nhl = normalize(ll - hl);
    // project into last cam basis
    vec2 luv = vec2(dot(nhl, r), dot(nhl, u));
    // project onto imaging plane NDC coords
    luv = luv / dot(nhl, lf) * FOV / vec2(asp, 1.0);
    
    if (any(greaterThan(luv, vec2(1.0))) || any(lessThan(luv, vec2(-1.0))))
        return vec4(0.0);
    
    // ndc to image pixel coords (minus half pixel, glsl uses 'pixel centers')
    vec2 fuv = (luv * -0.5 + 0.5) * iChanRes - 0.5;
    ivec2 iuv = ivec2(fuv);
    
    // samples with matching materials are considered
    vec4 col   = texelFetch(iChannel, iuv, 0);
    float c1 = float(int(decode(col.a) == ho));
    vec4 colx  = texelFetch(iChannel, iuv + ivec2(1, 0), 0);
    float c2 = float(int(decode(colx.a) == ho));
    vec4 coly  = texelFetch(iChannel, iuv + ivec2(0, 1), 0);
    float c3 = float(int(decode(coly.a) == ho));
    vec4 colxy = texelFetch(iChannel, iuv + ivec2(1, 1), 0);
    float c4 = float(int(decode(colxy.a) == ho));
    
    // fancy mixing
    vec2 duv = fuv - vec2(iuv);
    duv = duv*duv*(3.0 - 2.0*duv);
    return mix(mix(col*c1,colx*c2, duv.x), mix(coly*c3,colxy*c4, duv.x), duv.y);
}

vec4 render(in vec2 fragCoord, in vec2 iResolution, in int iFrame, in sampler2D lastFrame) {
    
    //lightLoc.x =  4.0 - cos(float(iFrame)*0.02);
    //lightLoc.y =  5.0 - sin(float(iFrame)*0.02);
    //lightLoc.z = -4.0 - sin(float(iFrame)*0.02);
    
    bool biased = false;
    float asp = iResolution.x / iResolution.y;
    vec2 ndc = (2.0 * fragCoord / iResolution - 1.0) * vec2(asp, 1.0);
    
    // get last cam
	vec3 ll = texelFetch(lastFrame, ivec2(0, int(iResolution.y) - 1), 0).xyz,
		 lo = texelFetch(lastFrame, ivec2(1, int(iResolution.y) - 1), 0).xyz,
	     ld = rotateXY(normalize(vec3(ndc, FOV)), lo.xy);
    
    // get current cam
	vec3 l = loc, // textureLod(controls, POS / iResolution, 0.0).xyz,
		 o = orient, // textureLod(controls, ROT / iResolution, 0.0).xyz,
         d = rotateXY(normalize(vec3(ndc, FOV)), o.xy),
         v = l - ll;
    
    // save current camera info
    if (int(fragCoord.y) == int(iResolution.y) - 1) {
        int x = int(fragCoord.x);
        if (x == 0) {
            return vec4(l, 0.0);
        } else if (x == 1) {
            return vec4(o, 0.0);
        }
    }
    
    // accumulated output
    vec3 outp = vec3(0.0);
    // current ray color
    vec3 rc = vec3(1.0);
    
    // reproject march
    vec3 res = march(l, d);
    vec3 hl = l + d * res.x;
    int ho = int(res.z);
    int _ho = ho;
    // reproject last frame onto this one
    vec4 fcol = reproject(ll, lo, hl, ho, asp, lastFrame, iResolution);
    if (fcol.a < eps) {
        // if sample failed, start with biased sampling to minimize variance
        biased = true;
    }
    
#ifdef BIASED
    // force biased rendering
    biased = true;
#endif
    
    if (fcol.a > float(TEMPORALSMOOTHING)) {
        fcol *= 1.0 - 1.0 / float(TEMPORALSMOOTHING);
    }
    
    for (int b = 0; b < BOUNCES; ++b) {
        if (b != 0) {
            res = march(l, d);
            ho = int(res.z);
            hl = l + d * res.x;
        }
        mat3 surf = getSurface(ho, hl);
        // emissive
        outp += rc * surf[1];
        // not nothing
        if (ho != LIGHT) {
            vec3 hn = norm(hl, eps);
            // diffuse color
            rc *= surf[0];
            vec3 smp = vec3(0.0);
            ivec3 s = seed(iFrame*iFrame, b, ivec2(fragCoord.xy*iResolution.yx));
            // do biased sampling
            if (biased) {
                // multiple importance sample lights and surfaces
                smp += MIS(hl, hn, ho, s);
                outp += rc * smp;
            } else {
                // unbiased sampling (GT)
                for (int i = 0; i < SMP_UNBIAS; ++i) {
                    vec3 sampleDir = cosHemiDir(hn, s + i);
                    float lpdf = cosHemiPDF(hn, sampleDir);
                    smp += lightMarch(hl, hn, sampleDir) / lpdf;
                }
                outp += rc * (smp / float(SMP_UNBIAS));
            }
            // modify ray for bounce
            brdf(l, d, hl, hn, s);
        }
    }
    outp /= float(BOUNCES);
    fcol.a = floor(fcol.a);
    fcol += vec4(outp, 1.0);
#ifndef BIASED
    if (biased) {
        fcol *= float(BIAS_WEIGHT);
    }
#endif
    fcol.a += encode(_ho);
    return fcol;
}

void main() {
    fragColor = render(gl_FragCoord.xy, iResolution, iFrame, iChannel0);
}
