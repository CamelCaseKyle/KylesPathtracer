////////////////////////////// Rendering Configuration ////////////////////////////// 

// always use biased sampling (fallback to unbiased ground truth)
#define BIASED
// ray bounces
#define BOUNCES 1
// raymarching steps
#define STEPS 255
// number of frames to reproject and smooth over time
#define TEMPORALSMOOTHING 16

// direct lit lambert surface
#define SMP_DIRECT_LAMBERT 1
// lambert surface lit lambert surface
#define SMP_LAMBERT_SURFACE_LAMBERT 1
// phong surface lit lambert surface
#define SMP_LAMBERT_SURFACE_PHONG 1

// direct lit phong surface
#define SMP_DIRECT_PHONG 1
// lambert surface lit phong surface
#define SMP_PHONG_SURFACE_LAMBERT 1
// phong surface lit phong surface
#define SMP_PHONG_SURFACE_PHONG 1

// unbiased light samples
#define SMP_UNBIAS 4
// bias sample weight to reduce initial unbiased variance
#define BIAS_WEIGHT 1

////////////////////////////// Math Toolkit ////////////////////////////// 

const float	eps = 0.001,     ieps = 0.999,   zfar = 50.0,       FOV = 1.5,
            HPI = 1.5707963, PI = 3.1415926, TWOPI = 6.2831853, SQRT2 = 1.4142136, SC45 = 0.7071068;
const vec2	  VEL = vec2(0.5, 0.5),   POS = vec2(1.5, 0.5),   ROT = vec2(2.5, 0.5),   MOU = vec2(3.5, 0.5);
const vec2	L_VEL = vec2(4.5, 0.5), L_POS = vec2(5.5, 0.5), L_ROT = vec2(6.5, 0.5), L_MOU = vec2(7.5, 0.5);

// generate a unique value for each pixel/frame
int genSeed(in int iFrame, in ivec2 fragCoord, in ivec2 iResolution) {
    return ( (iFrame<<12) + fragCoord.x + (fragCoord.y<<1) ) ^ fragCoord.x*iResolution.y ^ fragCoord.y*iResolution.x;
}

vec3 weyl3(in int v) {
    return fract(vec3(v*ivec3(13743434, 11258243, 9222443)) / exp2(24.0));
}

// fudged a bit to cut it off close to the (0,1) interval
vec3 logit3(in vec3 v) {
    vec3 t = 0.988 * (v + 0.006);
    return log(t / (1.0 - t)) * 0.221 + 0.5;
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

// thanks iq http:// iquilezles.org/www/articles/texture/texture.htm
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

float inRange(float a, float x, float b) {
    return step(a, x) * step(x, b);
}

// input is normalized visible spectrum wavelength (0. = 400nm violet/uv, 1. = 700nm red/ir)
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

// Thanks Paniq
vec3 linear_srgb(vec3 x) {
    return mix(1.055*pow(x, vec3(1./2.4)) - 0.055, 12.92*x, step(x, vec3(0.0031308)));
}

vec3 srgb_linear(vec3 x) {
    return mix(pow((x + 0.055)/1.055,vec3(2.4)), x / 12.92, step(x, vec3(0.04045)));
}

// Paniq's ACES fitted from https:// github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
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

float Lambertian(in vec3 hn, in vec3 nlv) {
    return max(eps, dot(nlv, hn));
}

float Phong(in vec3 rd, in vec3 hn, in vec3 nlv, in float gloss) {
    return pow(max(eps, dot(nlv, reflect(rd, hn))), gloss);
}

// using gaussian distribution
vec3 uniformSphere(int seed) {
    vec3 rnd = logit3(weyl3(seed));
    return rnd * 2.0 - 1.0;
}

// using normalized gaussian distribution
vec3 uniformDir(int seed) {
    return normalize(uniformSphere(seed));
}

// only one hemisphere
vec3 uniformHemiDir(vec3 hn, int seed) {
    vec3 rnd = uniformDir(seed);
    return rnd * sign(dot(hn, rnd));
}

// cosine distribution
vec3 cosHemiDir(vec3 hn, int seed) {
    vec3 rnd = uniformDir(seed);
    return normalize(hn + rnd * ieps);
}

// uniform sample cone
vec3 uniformConeDir(vec3 lv, float lr, in int seed) {
    vec3 rnd = weyl3(seed);
    float sa = linearAngle(length(lv), lr);
    float rad = sqrt(rnd.x) * tan(sa);
    float tha = rnd.y * TWOPI;
    vec3 r, u, nlv = normalize(lv);
    basis(nlv, r, u);
    return normalize(nlv + rad * (r * cos(tha) + u * sin(tha)));
}
////////////////////////////// Scene Modeling ////////////////////////////// 

vec2 sdMin(in vec2 a, in vec2 b) {
    if (a.x < b.x)
        return a;
    return b;
}

// Thanks iq and Dave Smith
float smin( float a, float b, float k ) {
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*k*0.25;
}

float smax(float a, float b, float k) {
    return -smin(-a,-b,k);
}

float sdBox(in vec3 p, in mat3 o, in vec3 s) {
    vec3 d = abs(p * o) - s;
    return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
}

const int
    LIGHT = 1,
    FLOOR = 2,
    WALL1 = 3,
    BOX = 4,
    WALL2 = 6,
    CEIL = 7;

// light parameters
vec4 light = vec4(6.0, 5.0, -4.0, 1.0);
vec3 lightColor = vec3(10.);
// plane parameters
vec4 flr = vec4( 0.0, 1.0, 0.0, 0.0);
vec4 cil = vec4( 0.0,-1.0, 0.0, 10.0);
vec4 wa1 = vec4(-1.0, 0.0, 0.0, 10.0);
vec4 wa2 = vec4( 0.0, 0.0, 1.0, 10.0);

mat3 getSurface(in int ho, in vec3 hl) {
    mat3 ret;
    if (ho == LIGHT) {
        ret[0] = vec3(1.0); // reflection color
        ret[1] = lightColor; // emission color
        ret[2] = vec3(1.0, 1.0, 0.0); // diffuse (x) specular (y)
    } else if (ho == BOX) {
        ret[0] = vec3(0.025 + 0.1 * float(int(floor(hl.x * 4.0) + floor(hl.y * 4.0) + floor(hl.z * 4.0)) & 1)); // reflection color
        ret[1] = vec3(0.0); // emission color
        ret[2] = vec3(1.0, 1.0, 0.0); // diffuse (x) specular (y)
    } else if (ho < 1) {
        
    } else {
        float refl = float(int(ho == FLOOR || ho == CEIL)) * (0.5 + float(int(floor(hl.x) + floor(hl.y) + floor(hl.z)) & 1)) * 0.2 + 0.8;
        float matCol = float(ho);
        float cm = cos(matCol) * 0.025;
        float sm = sin(matCol) * 0.025;
        ret[0] = vec3(0.05 + cm, 0.05 + sm, 0.05 - (cm + sm) * 0.25) * refl;
        ret[1] = vec3(0.0); // emission color
        ret[2] = vec3(refl, refl * 0.5, 0.0); // diffuse (x) specular (y)
    }
    // loose the average energy for the material, based on measurements from unbiased rendering
    // ret[2] *= 0.07;
    ret[2] *= 0.7;
    return ret;
}

vec2 sdf(in vec3 l, in int o) {
    vec2 d = vec2(zfar, 0.);
    if (o != FLOOR) d = sdMin(d, vec2(dot(l, flr.xyz) + flr.w, FLOOR));
    if (o != CEIL)  d = sdMin(d, vec2(dot(l, cil.xyz) + cil.w, CEIL));
    if (o != WALL1) d = sdMin(d, vec2(dot(l, wa1.xyz) + wa1.w, WALL1));
    if (o != WALL2) d = sdMin(d, vec2(dot(l, wa2.xyz) + wa2.w, WALL2));
    if (o != LIGHT) d = sdMin(d, vec2(length(l - light.xyz) - light.w, LIGHT));
    if (o != BOX)   d = sdMin(d, vec2(sdBox(l - vec3(7.5, 0.93, -7.5), mat3(1.0), vec3(0.8)) - 0.1, BOX));
    return d;
}

// Thanks Nimitz https:// www.shadertoy.com/view/Xts3WM
vec4 norcurv(in vec3 p, in float ep) {
    vec2 e = vec2(-1., 1.)*ep;
    float t1 = sdf(p + e.yxx, -1).x, t2 = sdf(p + e.xxy, -1).x;
    float t3 = sdf(p + e.xyx, -1).x, t4 = sdf(p + e.yyy, -1).x;
    return vec4(normalize(e.yxx*t1 + e.xxy*t2 + e.xyx*t3 + e.yyy*t4), .25/e.y*(t1 + t2 + t3 + t4 - 4.0*sdf(p, -1).x));
}

vec2 march(in vec3 l, in vec3 rd, in int o) {
    float t = 0.0;
    vec2 sdSmp;
    for (int i = 0; i < STEPS; ++i) {
        sdSmp = sdf(l + rd * t, o);
        t += sdSmp.x;
        if (sdSmp.x < eps)
            break;
        if (t > zfar)
            return vec2(zfar, 0.0);
    }
    return vec2(min(t, zfar), sdSmp.y);
}

////////////////////////////// PDF only ////////////////////////////// 

// surface illuminated by a sphere light, returns sample direction XYZ, pdf W
vec4 SphereLightPDF(in vec3 hl, in vec3 hn, in vec4 li, in int seed) {
    vec3 lv = li.xyz - hl;
    vec3 lambDir = uniformConeDir(lv, li.w, seed);
    float lpdf = solidAngle(dot(lv,lv), li.w*li.w);
    return vec4(lambDir, lpdf);
}

// surface illuminated by a lambertian surface, returns sample direction XYZ, pdf W
vec4 LambertPlanePDF(in vec3 hl, in vec3 hn, in vec4 li, in vec4 pl, in int seed) {
    // dot project light onto plane
    vec3 d = li.xyz - pl.xyz * (dot(li.xyz, pl.xyz) + pl.w);
    // surface to diffuse point d
    vec3 dv = d - hl;
    // diffuse point to light
    vec3 ld = li.xyz - d;
    // cosine sample disc on ground below
    float frad = min(length(dv), length(ld)) * 0.9;
    vec3 lambDir = uniformConeDir(dv, frad, seed);
    // PDF book keeping
    float lpdf = solidAngle(dot(dv,dv), frad*frad) / PI;
    float g2pdf = Lambertian(pl.xyz, -lambDir);
    return vec4(lambDir, lpdf * g2pdf);
}

// surface illuminated by a lambertian surface, returns sample direction XYZ, pdf W
vec4 PhongPlanePDF(in vec3 hl, in vec3 hn, in vec4 li, in vec4 pl, in int seed) {
    // get distance from each point to the plane
    float a = dot(hl, pl.xyz) + pl.w;
    float b = dot(li.xyz, pl.xyz) + pl.w;
    // similar triangles, just use the ratio of side lengths
    vec3 s = mix(hl - a*pl.xyz, li.xyz - b*pl.xyz, a / (a + b));
    // surface to specular point s
    vec3 sv = s - hl;
    float lsv = sqrt(dot(sv, sv)) * li.w;
    // specular point to light
    vec3 ls = li.xyz - s;
    // sample the image of the specular reflection on the surface TODO: implement roughness
    vec3 ts = sv * sqrt(dot(ls, ls));
    vec3 phongDir = uniformConeDir(ts, lsv, seed);
    // PDF book keeping
    float lpdf = solidAngle(dot(ts,ts), lsv*lsv) / PI;
    float spdf = Schlick(1.0, 3.0, dot(normalize(sv), pl.xyz));
    return vec4(phongDir, lpdf * spdf);
}

////////////////////////////// sampling only ////////////////////////////// 

// marches the light, returns contribution RGB
vec3 LightContribution(in vec3 hl, in int ho, in vec4 lvpdf) {
    vec2 lm = march(hl, lvpdf.xyz, ho);
    if (int(lm.y) == LIGHT)
        return lightColor * lvpdf.w;
    return vec3(0.0);
}

// marches the surface and light, returns contribution RGB
vec3 LambertPlaneContrib(in vec4 lvpdf, in vec3 hl, in int ho, in vec4 pl, in int po, in vec4 li, in int seed) {
    // diffuse plane sample
    vec2 dres = march(hl, lvpdf.xyz, ho);
    if (int(dres.y) != po)
        return vec3(0.0);
    // move off the surface a little bit
    vec3 _hl = hl + lvpdf.xyz * dres.x + pl.xyz * eps;
    vec3 _lv = li.xyz - _hl;
    float _lv2 = dot(_lv, _lv);
    // light sampling
    vec3 sampleDir = uniformConeDir(_lv, li.w, seed);
    vec3 lc = LightContribution(_hl, po, vec4(sampleDir, lvpdf.w));
    mat3 surf = getSurface(po, _hl);
    // surface emission, energy, color, light energy
    return surf[1] + surf[2].x * surf[0] * lc;
}

// marches the surface and light, returns contribution RGB
vec3 PhongPlaneContrib(in vec4 lvpdf, in vec3 hl, in int ho, in vec4 pl, in int po, in vec4 li, in int seed) {
    // specular plane sample
    vec2 sres = march(hl, lvpdf.xyz, ho);
    if (int(sres.y) != po)
        return vec3(0.0);
    // move off the surface a little bit
    vec3 _hl = hl + lvpdf.xyz * sres.x + pl.xyz * eps;
    vec3 _lv = li.xyz - _hl;
    float _lv2 = dot(_lv, _lv);
    // light sampling
    vec3 sampleDir = uniformConeDir(_lv, li.w, seed);
    vec3 lc = LightContribution(_hl, po, vec4(sampleDir, lvpdf.w));
    mat3 surf = getSurface(po, _hl);
    // surface emission, energy, color, light energy
    return surf[1] + surf[2].y * surf[0] * lc;
}

////////////////////////////// Sampling Strategy ////////////////////////////// 

// do unbiased sampling (diffuse BRDF)
vec3 UnbiasedLambertian(in vec3 hl, in vec3 hn, in int ho, in int seed) {
    // unbiased sampling (GT)
    vec3 smpUnbias = vec3(0.0);
    for (int i = 0; i < SMP_DIRECT_LAMBERT; ++i) {
        int si = seed + i;
        // sample light
        smpUnbias += LightContribution(hl, ho, vec4(cosHemiDir(hn, si), PI));
    }
    return smpUnbias /= float(SMP_DIRECT_LAMBERT);
}

// do unbiased sampling (specular BRDF)
vec3 UnbiasedPhong(in vec3 rd, in vec3 hl, in vec3 hn, in int ho, in int seed) {
    // unbiased sampling (GT)
    vec3 smpUnbias = vec3(0.0);
    for (int i = 0; i < SMP_DIRECT_LAMBERT; ++i) {
        int si = seed + i;
        // sample light
        smpUnbias += LightContribution(hl, ho, vec4(reflect(rd, hn), 1.0)); // not sure, replace with cone sample based on roughness
    }
    return smpUnbias /= float(SMP_DIRECT_LAMBERT);
}

// unbiased diffuse
void brdfLambertian(inout vec3 ro, inout vec3 rd, in vec3 hl, in vec3 hn, in int seed) {
    ro = hl + hn * eps;
    rd = cosHemiDir(hn, seed);
}

// unbiased specular
void brdfPhong(inout vec3 ro, inout vec3 rd, in vec3 hl, in vec3 hn, in ivec3 seed) {
    ro = hl + hn * eps;
    rd = reflect(rd, hn); // not sure, replace with cone sample based on roughness
}

// diffuse MIS
vec3 DMIS(in vec3 hl, in vec3 hn, in int ho, in int seed) {
    vec3 ret = vec3(0.0);

#if SMP_DIRECT_LAMBERT
    // take direct light samples
    vec3 smpDirect = vec3(0.0);
    for (int i = 0; i < SMP_DIRECT_LAMBERT; ++i) {
        int si = seed + i;
        // get all PDFs
        vec4 dlpdf = SphereLightPDF(hl, hn, light, si);
        // apply lambert pdf
        dlpdf.w *= Lambertian(hn, dlpdf.xyz);
        // compute CDF
        // float lightCDF = dlpdf.w;
        // roulette sampling
        // vec3 rnd = weyl3(si.z) * lightCDF;
        // if (rnd.z <= lightCDF) {
            smpDirect += LightContribution(hl, ho, dlpdf); // * (lightCDF / max(eps, dlpdf.w));
        // }
    }
    ret += smpDirect / float(SMP_DIRECT_LAMBERT);
#endif

#if SMP_LAMBERT_SURFACE_LAMBERT
    // take indirect diffuse samples
    vec3 smpIDD = vec3(0.0);
    for (int i = 0; i < SMP_LAMBERT_SURFACE_LAMBERT; ++i) {
        int si = seed + i;
        // get all PDFs
        vec4 fdpdf = LambertPlanePDF(hl, hn, light, flr, si);
        vec4 cdpdf = LambertPlanePDF(hl, hn, light, cil, si);
        vec4 w1dpdf = LambertPlanePDF(hl, hn, light, wa1, si);
        vec4 w2dpdf = LambertPlanePDF(hl, hn, light, wa2, si);
        // apply lambert pdf
        fdpdf.w *= Lambertian(hn, fdpdf.xyz);
        cdpdf.w *= Lambertian(hn, cdpdf.xyz);
        w1dpdf.w *= Lambertian(hn, w1dpdf.xyz);
        w2dpdf.w *= Lambertian(hn, w2dpdf.xyz);
        // compute CDF
        float floorDiffCDF = fdpdf.w;
        float ceilDiffCDF = floorDiffCDF + cdpdf.w;
        float w1DiffCDF = ceilDiffCDF + w1dpdf.w;
        float w2DiffCDF = w1DiffCDF + w2dpdf.w;
        // roulette sampling
        vec3 rnd = weyl3(si) * w2DiffCDF;
        if (rnd.z <= floorDiffCDF)
            smpIDD += LambertPlaneContrib(fdpdf, hl, ho, flr, FLOOR, light, si) * (w2DiffCDF / max(eps, fdpdf.w));
        else if (rnd.z <= ceilDiffCDF)
            smpIDD += LambertPlaneContrib(cdpdf, hl, ho, cil, CEIL, light, si) * (w2DiffCDF / max(eps, cdpdf.w));
        else if (rnd.z <= w1DiffCDF)
            smpIDD += LambertPlaneContrib(w1dpdf, hl, ho, wa1, WALL1, light, si) * (w2DiffCDF / max(eps, w1dpdf.w));
        else
            smpIDD += LambertPlaneContrib(w2dpdf, hl, ho, wa2, WALL2, light, si) * (w2DiffCDF / max(eps, w2dpdf.w));
    }
    ret += smpIDD / float(SMP_LAMBERT_SURFACE_LAMBERT);
#endif

#if SMP_LAMBERT_SURFACE_PHONG
    // take indirect specular samples
    vec3 smpIDS = vec3(0.0);
    for (int i = 0; i < SMP_LAMBERT_SURFACE_PHONG; ++i) {
        int si = seed + i;
        // get all PDFs
        vec4 fspdf = PhongPlanePDF(hl, hn, light, flr, si);
        vec4 cspdf = PhongPlanePDF(hl, hn, light, cil, si);
        vec4 w1spdf = PhongPlanePDF(hl, hn, light, wa1, si);
        vec4 w2spdf = PhongPlanePDF(hl, hn, light, wa2, si);
        // apply lambert pdf
        fspdf.w *= Lambertian(hn, fspdf.xyz);
        cspdf.w *= Lambertian(hn, cspdf.xyz);
        w1spdf.w *= Lambertian(hn, w1spdf.xyz);
        w2spdf.w *= Lambertian(hn, w2spdf.xyz);
        // compute CDF
        float floorSpecCDF = fspdf.w;
        float ceilSpecCDF = floorSpecCDF + cspdf.w;
        float w1SpecCDF = ceilSpecCDF + w1spdf.w;
        float w2SpecCDF = w1SpecCDF + w2spdf.w;
        // roulette sampling
        vec3 rnd = weyl3(si) * w2SpecCDF;
        if (rnd.z <= floorSpecCDF)
            smpIDS += PhongPlaneContrib(fspdf, hl, ho, flr, FLOOR, light, si) * (w2SpecCDF / max(eps, fspdf.w));
        else if (rnd.z <= ceilSpecCDF)
            smpIDS += PhongPlaneContrib(cspdf, hl, ho, cil, CEIL, light, si) * (w2SpecCDF / max(eps, cspdf.w));
        else if (rnd.z <= w1SpecCDF)
            smpIDS += PhongPlaneContrib(w1spdf, hl, ho, wa1, WALL1, light, si) * (w2SpecCDF / max(eps, w1spdf.w));
        else
            smpIDS += PhongPlaneContrib(w2spdf, hl, ho, wa2, WALL2, light, si) * (w2SpecCDF / max(eps, w2spdf.w));
    }
    ret += smpIDS / float(SMP_LAMBERT_SURFACE_PHONG);
#endif

    return ret;
}

// specular MIS
vec3 SMIS(in vec3 rd, in vec3 hl, in vec3 hn, in int ho, in int seed) {
    vec3 ret = vec3(0.0);

#if SMP_DIRECT_PHONG
    // take direct light samples
    vec3 smpDirect = vec3(0.0);
    for (int i = 0; i < SMP_DIRECT_PHONG; ++i) {
        int si = seed + i;
        // get all PDFs
        vec4 dlpdf = SphereLightPDF(hl, hn, light, si);
        // apply phong pdf
        dlpdf.w *= Phong(rd, hn, dlpdf.xyz, 5.0);
        // compute CDF
        // float lightCDF = dlpdf.w;
        // roulette sampling
        // vec3 rnd = weyl3(si) * lightCDF;
        // if (rnd.z <= lightCDF) {
            smpDirect += LightContribution(hl, ho, dlpdf); // * (lightCDF / max(eps, dlpdf.p));
        // }
    }
    ret += smpDirect / float(SMP_DIRECT_PHONG);
#endif

#if SMP_PHONG_SURFACE_LAMBERT
    // take indirect diffuse samples
    vec3 smpIDD = vec3(0.0);
    for (int i = 0; i < SMP_PHONG_SURFACE_LAMBERT; ++i) {
        int si = seed + i;
        // get all PDFs
        vec4 fdpdf = LambertPlanePDF(hl, hn, light, flr, si);
        vec4 cdpdf = LambertPlanePDF(hl, hn, light, cil, si);
        vec4 w1dpdf = LambertPlanePDF(hl, hn, light, wa1, si);
        vec4 w2dpdf = LambertPlanePDF(hl, hn, light, wa2, si);
        // apply phong pdf
        fdpdf.w *= Phong(rd, hn, fdpdf.xyz, 5.0);
        cdpdf.w *= Phong(rd, hn, cdpdf.xyz, 5.0);
        w1dpdf.w *= Phong(rd, hn, w1dpdf.xyz, 5.0);
        w2dpdf.w *= Phong(rd, hn, w2dpdf.xyz, 5.0);
        // compute CDF
        float floorDiffCDF = fdpdf.w;
        float ceilDiffCDF = floorDiffCDF + cdpdf.w;
        float w1DiffCDF = ceilDiffCDF + w1dpdf.w;
        float w2DiffCDF = w1DiffCDF + w2dpdf.w;
        // roulette sampling
        vec3 rnd = weyl3(si) * w2DiffCDF;
        if (rnd.z <= floorDiffCDF)
            smpIDD += LambertPlaneContrib(fdpdf, hl, ho, flr, FLOOR, light, si) * (w2DiffCDF / max(eps, fdpdf.w));
        else if (rnd.z <= ceilDiffCDF)
            smpIDD += LambertPlaneContrib(cdpdf, hl, ho, cil, CEIL, light, si) * (w2DiffCDF / max(eps, cdpdf.w));
        else if (rnd.z <= w1DiffCDF)
            smpIDD += LambertPlaneContrib(w1dpdf, hl, ho, wa1, WALL1, light, si) * (w2DiffCDF / max(eps, w1dpdf.w));
        else
            smpIDD += LambertPlaneContrib(w2dpdf, hl, ho, wa2, WALL2, light, si) * (w2DiffCDF / max(eps, w2dpdf.w));
    }
    ret += smpIDD / float(SMP_PHONG_SURFACE_LAMBERT);
#endif

#if SMP_PHONG_SURFACE_PHONG
    // take indirect specular samples
    vec3 smpIDS = vec3(0.0);
    for (int i = 0; i < SMP_PHONG_SURFACE_PHONG; ++i) {
        int si = seed + i;
        // get all PDFs
        vec4 fspdf = PhongPlanePDF(hl, hn, light, flr, si);
        vec4 cspdf = PhongPlanePDF(hl, hn, light, cil, si);
        vec4 w1spdf = PhongPlanePDF(hl, hn, light, wa1, si);
        vec4 w2spdf = PhongPlanePDF(hl, hn, light, wa2, si);
        // apply phong pdf
        fspdf.w *= Phong(rd, hn, fspdf.xyz, 5.0);
        cspdf.w *= Phong(rd, hn, cspdf.xyz, 5.0);
        w1spdf.w *= Phong(rd, hn, w1spdf.xyz, 5.0);
        w2spdf.w *= Phong(rd, hn, w2spdf.xyz, 5.0);
        // compute CDF
        float floorSpecCDF = fspdf.w;
        float ceilSpecCDF = floorSpecCDF + cspdf.w;
        float w1SpecCDF = ceilSpecCDF + w1spdf.w;
        float w2SpecCDF = w1SpecCDF + w2spdf.w;
        vec3 rnd = weyl3(si) * w2SpecCDF;
        if (rnd.z <= floorSpecCDF)
            smpIDS += PhongPlaneContrib(fspdf, hl, ho, flr, FLOOR, light, si) * (w2SpecCDF / max(eps, fspdf.w));
        else if (rnd.z <= ceilSpecCDF)
            smpIDS += PhongPlaneContrib(cspdf, hl, ho, cil, CEIL, light, si) * (w2SpecCDF / max(eps, cspdf.w));
        else if (rnd.z <= w1SpecCDF)
            smpIDS += PhongPlaneContrib(w1spdf, hl, ho, wa1, WALL1, light, si) * (w2SpecCDF / max(eps, w1spdf.w));
        else
            smpIDS += PhongPlaneContrib(w2spdf, hl, ho, wa2, WALL2, light, si) * (w2SpecCDF / max(eps, w2spdf.w));
    }
    ret += smpIDS / float(SMP_PHONG_SURFACE_PHONG);
#endif

    return ret;
}
////////////////////////////// Buffers ////////////////////////////// 

void encodeGbuffer(out vec4 val, in vec3 n, in int o, in float t) {
    val = vec4(n * float(o), t);
}

void decodeGbuffer(in vec4 val, out vec3 n, out int o, out float t) {
    t = val.w;
    o = int(length(val.xyz)+eps);
    n = val.xyz / float(o);
}

float encodeBuffer(in int o) {
    return float(o) * 0.01001;
}

int decodeBuffer(in float a) {
    return int(fract(a) * 100.1);
}

// use decodeAll macro below instead
void _decodeAll(in vec2 fragCoord, in vec2 iResolution, in int iFrame, in sampler2D gBuffer, in sampler2D lastFrame, out float asp, out vec3 rl, out vec2 ro, out vec3 rd, out vec3 ll, out vec2 lo, out vec3 ld, out float vv, out vec3 hn, out int ho, out vec3 hl) {
    asp = iResolution.x / iResolution.y;
    vec2 ndca = (2.0 * fragCoord.xy / iResolution.xy - 1.0) * vec2(asp, 1.0);
    // get current cam
    if (iFrame > 1) {
        rl = texelFetch(gBuffer, ivec2(0, int(iResolution.y) - 1), 0).xyz;
        ro = texelFetch(gBuffer, ivec2(1, int(iResolution.y) - 1), 0).xy;
        ll = texelFetch(lastFrame, ivec2(0, int(iResolution.y) - 1), 0).xyz;
        lo = texelFetch(lastFrame, ivec2(1, int(iResolution.y) - 1), 0).xy;
    }
    rd = rotateXY(normalize(vec3(ndca, FOV)), ro);
    ld = rotateXY(normalize(vec3(ndca, FOV)), lo);
	vv = length(rl - ll);
    // get last trace
    float ht;
    decodeGbuffer(texelFetch(gBuffer, ivec2(fragCoord.xy), 0), hn, ho, ht);
    hl = rl + rd * ht;
}

// defines variables for current ray, last ray, current hit, velocity, etc
#define decodeAll(gBuffer, lastFrame) int ho; float asp, vv; vec2 ro = orient.xy, lo = orient.xy; vec3 rl = loc.xyz, ll = loc.xyz, rd, ld, hn, hl; _decodeAll(fragCoord.xy, iResolution.xy, iFrame, gBuffer, lastFrame, asp, rl, ro, rd, ll, lo, ld, vv, hn, ho, hl)

// last loc, last orient, hit loc, hit object, 1/aspect, channel, resolution
vec4 reproject(in vec3 ll, in vec2 lo, in vec3 hl, in int ho, in float asp, sampler2D iChannel, in vec2 iChanRes) {
    // last camera basis
    vec3 lf = rotateXY(vec3(0.0, 0.0, 1.0), lo);
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
    float c1 = float(int(decodeBuffer(col.a) == ho));
    vec4 colx  = texelFetch(iChannel, iuv + ivec2(1, 0), 0);
    float c2 = float(int(decodeBuffer(colx.a) == ho));
    vec4 coly  = texelFetch(iChannel, iuv + ivec2(0, 1), 0);
    float c3 = float(int(decodeBuffer(coly.a) == ho));
    vec4 colxy = texelFetch(iChannel, iuv + ivec2(1, 1), 0);
    float c4 = float(int(decodeBuffer(colxy.a) == ho));
    
    // fancy mixing
    vec2 duv = fuv - vec2(iuv);
    // duv = duv*duv*(3.0 - 2.0*duv);
    return mix(mix(col*c1,colx*c2, duv.x), mix(coly*c3,colxy*c4, duv.x), duv.y);
}
