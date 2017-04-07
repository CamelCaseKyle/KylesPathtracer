#version 330 core

layout(location = 0) out vec4 fragColor;

in vec2 uv;

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

//#define SHOWSEG
#define MULTI		8
#define BOUNCE_PXL	6
#define GAMMA		.5
#define EXPOSURE	1.

#define tex(a,b) textureLod(a,b,0.)

#define vec3max(a) max(a.x, max(a.y, a.z))
#define vec2max(a) max(a.x, a.y)
#define vec3min(a) min(a.x, min(a.y, a.z))
#define vec2min(a) min(a.x, a.y)

#define valid(a) (a.t.y < zfar)
#define inRange(a,x,b) (step(a, x) * step(x, b))
#define overlap(a, b) ((a.t.y > b.t.x) && (a.t.x < b.t.y))
#define contains(a, b) ((a.t.x > b.t.x) && (a.t.y < b.t.y))
#define minT(a) (a.y<a.x)? zfar: (a.x<0.)? (a.y<0.)? zfar: a.y: a.x

const float zfar = 1000.,
			eps  = .00001, ieps = .99999,
			sml  = .001, isml = .999,
			sc45 = .7071067,
			pi_5 = 1.5707963, pi2  = 6.2831853,
			pi   = 3.1415926, pi_rcp = .3183098;

//Absorb, Emit rgb wavelengths, Surface scatter, sUbsurface scatter, Refractive index, Emission Uniformity, unique iDentifier
struct mat { vec3 a, e; vec2 s, u; float r, eu; int d; };
//1D line segment with ids for each end point
struct seg { vec2 t; ivec2 d; };
//Location, Normal, Info(distance, thickness, nothing), material
struct hit { vec3 l, n; float d; mat m; };
//Origin, Direction, Color, current material
struct ray { vec3 o, d, c; mat m; };
//center Location, Radius, Material
struct sph { vec3 l; float r; mat m; int d; };
//Location, Normal, Info(size x,y,'depth'), Material
struct pln { vec3 l, n, i; mat m; int d; };
//center, size, Material
struct box { vec3 c, s; vec2 t; mat m; int d; };

const mat nullMat = mat(vec3(0.), vec3(0.), vec2(0.), vec2(0.), 0., 0., 0);
const seg nullSeg = seg(vec2(zfar), ivec2(0));
const hit nullHit = hit(vec3(0.), vec3(0.), zfar, nullMat);
const vec2 nullT = vec2(zfar, 0.);

////////////////////// Materials
const mat
//conductors
    _air = mat(vec3(.99), vec3(0.), vec2(0.), vec2(0.1), 1.0003, 0., 0),
	_dmn = mat(vec3(.99), vec3(0.), vec2(0.), vec2(0.), 2.45, 0., 1),
//reflective insulator
    _chr = mat(vec3(.99), vec3(0.), vec2(0.), vec2(0.), -1., 0., 10),
    _gld = mat(vec3(.99, .9, 0.), vec3(0.), vec2(0.), vec2(0.), -1., 0., 11),
//diffuse insulator
    _wht = mat(vec3(.9), vec3(0.), vec2(1.), vec2(.2), -1., 0., 20),
    _cyn = mat(vec3(.1,.8,.9), vec3(0.), vec2(1.), vec2(.2), -1., 0., 21),
    _orn = mat(vec3(.9,.2,.1), vec3(0.), vec2(1.), vec2(.2), -1., 0., 22),
	_blk = mat(vec3(.01), vec3(0.), vec2(1.), vec2(.2), -1., 0., 23),
//emissive
	_wht_e = mat(vec3(.99), vec3(99.), vec2(1.), vec2(.2), -1., 0., 30),
//dynamic materials
	_mta = mat(vec3(0.), vec3(0.), vec2(0.), vec2(0.), 0., 0., 40),
    _nrm = mat(vec3(0.), vec3(0.), vec2(0.), vec2(0.), 0., 0., 41);

////////////////////// Primitives
sph sph0 = sph(vec3(0., 0., 0.), 100., _cyn, 1),
//lights
	lit0 = sph(vec3(0., 25., 25.), 5., _wht_e, 9);
//planes
pln pln0 = pln(vec3(0., -5., 0.), vec3(0., 1., 0.), vec3(1., 1., .1), _wht, 10);
//lens shape
sph sph1 = sph(vec3(0.,0.,0.), 20., _dmn, 2),
	sph3 = sph(vec3(0.,25.,25.), 20., _dmn, 4);
//room shape
box box1 = box(vec3(0.), vec3(5.), vec2(0.), _blk, 21),
	box2 = box(vec3(0.), vec3(4.5), vec2(0.), _blk, 22),
	box3 = box(vec3(0.,3.,4.75), vec3(.71,1.,1.), vec2(0.), _blk, 23);
sph sph2 = sph(vec3(0., -28., -28.), 4., _gld, 3);
box box4 = box(vec3(0., -29., -29.), vec3(6.,.1,6.), vec2(0.), _chr, 24);

////////////////////// random number (distribution sucks)
float hash12( in vec2 p ) {
	float h = dot(p,vec2(101.7, 683.11));
    return fract(sin(h)*467.709);
}
float hash13( in vec3 p ) {
	float h = dot(p,vec3(571.127, 467.311, 881.521));
    return fract(sin(h)*467.197);
}
vec2 hash22( in vec2 p ) {
	vec2 h = vec2(dot(p,vec2(467.127, 881.311)),
                  dot(p,vec2(7.101, 11.683)));
    return fract(sin(h)*467.281);
}
vec3 hash33( in vec3 p ) {
	vec3 h = vec3(dot(p,vec3(467.127, 881.311, 571.521)),
                  dot(p,vec3(7.101, 11.683, 13.331)),
                  dot(p,vec3(683.331, 761.101, 823.127)));
    return fract(sin(h)*467.281);
}

////////////////////// normal(ish) distribution
float bell12( in vec2 n ) {
	float r0 = hash12(n + 7.07),
	 r1 = hash12(n + 11.11),
     r2 = hash12(n + 17.17);
	return (r0+r1+r2) * .3333;
}
float bell13( in vec3 n ) {
	float r0 = hash13(n + 7.07),
	 r1 = hash13(n + 11.11),
     r2 = hash13(n + 17.17);
	return (r0+r1+r2) * .3333;
}
vec3 bell33( in vec3 n ) {
	vec3 r0 = hash33(n + 7.07),
	 r1 = hash33(n + 11.11),
     r2 = hash33(n + 17.17);
	return (r0+r1+r2) * .3333;
}

////////////////////// perlin noise with scale, smoothsteped to 0. - 1.
float fade(float t) {
    return t * t * t * (t * (t * 6. - 15.) + 10.);
}
float pnoise12(vec2 p, float scl) {
    vec2 i = floor(p*scl),
     f = fract(p*scl),
     u = vec2(fade(f.x), fade(f.y)),
     o = vec2(0., 1.),
     g00 = hash22(i).xy, g01 = hash22(i + o.xy).xy,
     g11 = hash22(i + o.yy).xy, g10 = hash22(i + o.yx).xy,
     d00 = f, d01 = f - o.xy,
     d11 = f - o.yy, d10 = f - o.yx;
    float s00 = dot(g00, d00), s01 = dot(g01, d01),
     s11 = dot(g11, d11), s10 = dot(g10, d10),
     x1 = mix(s01,s11, u.x), x2 = mix(s00,s10, u.x);
    return mix(x2, x1, u.y);
}

////////////////////// some tools
vec3 rotateY(vec3 p, float angle) {
    float c = cos(angle), s = sin(angle);
    return vec3(c*p.x + s*p.z, p.y, -s*p.x + c*p.z);
}
vec3 rotateXY(in vec3 p, in vec2 angle) {
    vec2 c=cos(angle), s=sin(angle); vec3 o = p;
    o.yz *= mat2(c.x,s.x,-s.x,c.x);  o.xz *= mat2(c.y,s.y,-s.y,c.y);
    return o;
}
void basis(in vec3 n, out vec3 f, out vec3 r) {
    if(n.z < -0.999999) {
        f = vec3(0.,-1.,0.);
        r = vec3(-1.,0.,0.);
    } else {
    	float a = 1./(1. + n.z);
    	float b = -n.x*n.y*a;
    	f = vec3(1. - n.x*n.x*a, b, -n.x);
    	r = vec3(b, 1. - n.y*n.y*a , -n.y);
    }
}
vec3 slerp(in vec3 start, in vec3 end, in float percent) {
    float dt = dot(start, end), theta = acos(dt)*percent;
    return start*cos(theta) + normalize(end - start*dt)*sin(theta);
}
vec2 rap(in vec3 n) {
    return vec2(atan(n.z, n.x) + pi, acos(-n.y)) / vec2(pi2, pi);
}
vec3 rndSph(in vec3 chaos) {
    return normalize(bell33(chaos)*2.-1.);
}
vec3 rndHemi(in vec3 norm, in vec3 chaos) {
    vec3 guess = rndSph(chaos);
    return guess*(step(0., dot(norm, guess)) * 2. - 1.);
}
float schlick(in float r1, in float r2, in float vn) {
    float r0 = (r1 - r2) / (r1 + r2);
	return mix(r0*r0, 1., pow(1. - vn, 5.));
}
float getRoughness(in vec3 rd, in mat m, in float rnd) {
    return clamp(abs(bell13(rd + rnd) * 2. - 1.) * m.s.y + m.s.x, 0., 1.);
}
float getEmission(in vec3 rd, in vec3 nrm, in mat m) {
	return pow(clamp(-dot(nrm, rd), m.eu, 1.), 2.);   
}

//////////////////////Dynamic material implementation
mat mta(inout hit h, in ray r) {
    hash13(r.d + h.n + iGlobalTime);
    //use roughness
	h.n = -r.d;
    //retro reflective
    return _mta;
}
mat nrm(inout hit h, in ray r) {
    vec3 maxnorm = max(vec3(0.),h.n);
    return mat(vec3(.01), maxnorm*.1, vec2(1.), vec2(.5), -1., 0., _nrm.d); 
}

//gets dynamic material properties
void updateMaterials(inout hit ret, in ray r) {
    if (ret.m.d < 40) return;
    else if (ret.m.d == _mta.d) ret.m = mta(ret, r);
    else if (ret.m.d == _nrm.d) ret.m = nrm(ret, r);
}

////////////////////// SD functions
float sdSphere(in vec3 l, in sph s) {
    vec3 oc = l - s.l;
    return dot(oc, oc) - s.r * s.r;
}
float sdPln(in vec3 l, in pln p) {
    return dot(p.n, p.l - l);
}
float sdBox(in vec3 l, in box b) {
    vec3 d = abs(b.c - l) - b.s;
    return min(vec3max(d), 0.) + length(max(d, 0.));
}

////////////////////// Normal functions
vec3 nSphere(in vec3 l, in sph s) {
    return (l - s.l) / s.r;   
}
vec3 nPlane(in pln p) {
    return p.n;
}
vec3 nBox(in vec3 l, in box b) {
    vec3 a = l - b.c;
    return step(b.s*ieps, abs(a)) * sign(a);
}

////////////////////// Segment operators
vec2 lt(in seg s) {
    if (s.t.x < s.t.y && s.t.x > 0.) return vec2(s.t.x, float(s.d.x));
    else if (s.t.y > 0.) return vec2(s.t.y, float(s.d.y));
    return nullT;
}
vec2 pickXLT(in seg l, in seg r) {
    if (l.t.x < r.t.x) return vec2(l.t.x, float(l.d.x));
    return vec2(r.t.x, float(r.d.x));
}
vec2 pickYLT(in seg l, in seg r) {
    if (l.t.y < r.t.y) return vec2(l.t.y, float(l.d.y));
    return vec2(r.t.y, float(r.d.y));
}
vec2 pickXGT(in seg l, in seg r) {
    if (l.t.x > r.t.x && l.t.x < zfar) return vec2(l.t.x, float(l.d.x));
    else if (r.t.x < zfar) return vec2(r.t.x, float(r.d.x));
    return nullT;
}
vec2 pickYGT(in seg l, in seg r) {
    if (l.t.y > r.t.y && l.t.y < zfar) return vec2(l.t.y, float(l.d.y));
    else if (r.t.y < zfar) return vec2(r.t.y, float(r.d.y));
    return nullT;
}
void tIntersect(in seg l, in seg r, out seg o) {
    o = nullSeg;
	if (!overlap(l, r) || !valid(l) || !valid(r)) return;
    vec2 t1 = pickXGT(l, r), t2 = pickYLT(l, r);
    o = seg(vec2(t1.x, t2.x), ivec2(int(t1.y), int(t2.y)));
}
void tUnion(in seg l, in seg r, out seg o, out seg p) {
	o = nullSeg; p = nullSeg;
    if (!overlap(l, r)) { o = l; p = r; return; }
    vec2 t1 = pickXLT(l, r), t2 = pickYGT(l, r);
    o = seg(vec2(t1.x, t2.x), ivec2(int(t1.y), int(t2.y)));
}
void tDiff(in seg l, in seg r, out seg o, out seg p) {
    o = nullSeg; p = nullSeg;
    if (!overlap(l, r)) { o = l; return; }
    if (contains(l, r)) return;
    if (contains(r, l)) {
		vec2 t1 = pickXLT(l, r), t2 = pickXGT(l, r), t3 = pickYLT(l, r), t4 = pickYGT(l, r);
	    o = seg(vec2(t1.x, t2.x), ivec2(int(t1.y), int(t2.y)));
	    p = seg(vec2(t3.x, t4.x), ivec2(int(t3.y), int(t4.y)));
        return;
    }
    vec2 t1 = nullT, t2 = nullT;
    if (l.t.x < r.t.x) { t1 = pickXLT(l, r); t2 = pickXGT(l, r);
    } else { t1 = pickYLT(l, r); t2 = pickYGT(l, r); } 
    o = seg(vec2(t1.x, t2.x), ivec2(int(t1.y), int(t2.y)));
}

////////////////////// Segment functions
seg lt(in seg l, in seg r) {
    if (l.t.x < r.t.x && l.t.y > 0.) return l;
    else if (r.t.y > 0.) return r;
    return nullSeg;
}
seg tSphere(in ray r, in sph s) {
    vec3 oc = r.o - s.l;
    float c = dot(oc, oc) - s.r * s.r,
          b = -dot(oc, r.d),
          h = b*b - c;
    if (h < 0.) return nullSeg;
    h = sqrt(h);
    return seg(vec2(b-h,b+h), ivec2(s.d, -s.d));
}
seg tPlane(in ray r, in pln p) {
    float t = dot(p.n, p.l - r.o) / dot(p.n, r.d);
    return seg(vec2(t, t+eps), ivec2(p.d, -p.d));
}
seg tBox(in ray r, in box b) {
    vec3 t1 = (b.c-b.s - r.o)/r.d,
         t2 = (b.c+b.s - r.o)/r.d;
    float tn = vec3max(min(t1, t2)),
          tx = vec3min(max(t1, t2));
    if (tx<tn || tx<0.) return nullSeg;
    return seg(vec2(tn, tx), ivec2(b.d, -b.d));
}

////////////////////// Hit functions
hit lt(in hit l, in hit r) {
    if (l.d < r.d && l.d > 0.) return l;
    return r;
}
hit traceSphere(in ray r, in sph s) {
    seg g = tSphere(r, s);
    float d = minT(g.t);
    vec3 l = r.o + r.d * d;
    return hit(l, nSphere(l, s), d, s.m);   
}
hit tracePlane(in ray r, in pln p) {
    seg s = tPlane(r, p);
    float d = minT(s.t);
    return hit(r.o + r.d * d, p.n, d, p.m);
}
hit traceBox(in ray r, in box b) {
    seg s = tBox(r, b);
    float d = minT(s.t);
    vec3 l = r.o + r.d * d;
    return hit(l, nBox(l, b), d, b.m);
}
hit traceBoxTransformed(in ray r, in box b) {
    ray _r = ray(rotateXY(r.o - b.c, b.t), rotateXY(r.d, b.t), r.c, r.m);
    seg s = tBox(_r, b);
    float d = minT(s.t);
    vec3 l = _r.o + _r.d * d;
    return hit(l, rotateXY(nBox(l, b), -b.t), d, b.m);
}
hit lens1(in ray r) {
	//get segments of shapes
    seg a = tSphere(r, sph1), b = tSphere(r, sph3), c = nullSeg, d = nullSeg;
    tIntersect(a,b, c);
	sph1.l -= vec3(0.,6.,6.);
	sph3.l -= vec3(0.,6.,6.);
	a = tSphere(r, sph1); b = tSphere(r, sph3);
    tIntersect(a,b, d);
	sph1.l += vec3(0.,6.,6.);
	sph3.l += vec3(0.,6.,6.);
    //choose lowest non negitive point
    vec2 g = lt(lt(c,d));
    //if t point object is not defined
    if (g.y == 0.) return nullHit;
    //point properties
    mat o = nullMat;
    g.x = minT(vec2(g.x, zfar));
    vec3 n = vec3(0.), l = r.o + r.d * g.x;
    //get id of chosen hit and 'first' hit
    int id = abs(int(g.y));
    //choose material
    if (id == sph1.d) {
        o = sph1.m; n = nSphere(l, sph1);
    } else if (id == sph3.d) {
        o = sph3.m; n = nSphere(l, sph3);
	} else return nullHit;
    //invert normal if t is the second intersection
    return hit(l, n * sign(g.y), g.x, o);
}
hit room1(in ray r) {
	//get segments of shapes
    seg a = tBox(r, box1), b = tBox(r, box2), c = tBox(r, box3);
    seg[] s = seg[6](nullSeg, nullSeg, nullSeg, nullSeg, nullSeg, nullSeg);
    tDiff(a,b, s[0], s[1]);
    tDiff(s[0],c, s[2], s[3]);
    tDiff(s[1],c, s[4], s[5]);
    //choose lowest non negitive point
    vec2 g = lt(lt(s[2], lt(s[3], lt(s[4], s[5]))));
    //if t point object is not defined
    if (g.y == 0.) return nullHit;
    //point properties
    mat o = nullMat;
    g.x = minT(vec2(g.x, zfar));
    vec3 n = vec3(0.), l = r.o + r.d * g.x;
    //get id of chosen hit and 'first' hit
    int id = abs(int(g.y));
    //choose material
    if (id == box1.d) {
        o = box1.m; n = nBox(l, box1);
    } else if (id == box2.d) {
        o = box2.m; n = nBox(l, box2);
    } else if (id == box3.d) {
        o = box3.m; n = nBox(l, box3);
	} else return nullHit;
    //invert normal if t is the second intersection
    return hit(l, n * sign(g.y), g.x, o);
}

//calculates where the ray intersects
hit traceScene(in ray r) {
    hit ret = nullHit;
    
    ret = lt(ret, traceBox(r, box4));
    ret = lt(ret, traceSphere(r, sph1));
    ret = lt(ret, traceSphere(r, lit0));
    ret = lt(ret, traceSphere(r, sph0));
    //ret = lt(ret, traceSphere(r, sph2));
	
    ret = lt(ret, room1(r));
    //ret = lt(ret, lens1(r));
    
    //mirror normal if needed
    float nmod = step(dot(r.d, ret.n), 0.) * 2. - 1.;
	ret.n *= nmod;
    ret.l += ret.n * eps;
	
    return ret;
}

//calculates if the camera is inside anything for volumetrics
mat mapScene(in vec3 l) {
    hit h = nullHit;
    h = lt(h, hit(l, l, sdSphere(l, sph0), sph0.m));
    //h = lt(h, hit(l, l, vec3(sdSphere(l, sph1)), sph1.m));
	return h.m;
}

//ray reflect scatter
void brdf(inout ray r, inout hit res, in float rnd) {
    //get random point between ro and dist
    float d = hash13(r.d - rnd)*res.d,
        //how much the material changes the ray
        roughness = getRoughness(r.d, r.m, rnd);
    //new origin is randomly along ray
    r.o = r.o + r.d*d;
    //randomly reflect anywhere
    r.d = rndSph(res.l + rnd);
    r.c *= r.m.a;
}

//ray transmit scatter
void btdf(inout ray r, inout hit res, in float rnd) {
    //get random point between ro and dist
    float d = hash13(r.d - rnd)*res.d,
        //how much the air changes the ray
        roughness = getRoughness(r.d, r.m, rnd);
    //new origin is randomly along ray
    r.o = r.o + r.d*d;
    //new direction is randomly is cos weighted about r.d
    r.d = slerp(r.d, rndHemi(r.d , res.l + rnd), roughness);
    //r.d = slerp(r.d, rndHemi(r.d, r.o - rnd), roughness);
    r.c *= r.m.a;
}

//called at the interface of two materials
void bsdf(inout ray r, inout hit res, in float rnd) {
    //calculate roughness (randomness) of surface
    float roughness = getRoughness(r.d, res.m, rnd),
		  //decide the ray's fate
          sch = schlick(r.m.r, res.m.r, dot(res.n, -r.d)),
		  //reflect [1] or refract [-1]
		  ror = step(hash13(r.d - rnd), sch) * 2. - 1.;
	vec3 newDir = (ror < 0.)? refract(r.d, res.n, r.m.r / res.m.r): reflect(r.d, res.n);
    //new ray properties (advance the ray to the proper side of the surface)
	r.o = res.l + res.n*eps*ror;
    //new hemi about res.n weignted towards newDir
    r.d = slerp(newDir, rndHemi(res.n , res.l + rnd), roughness);
    r.c *= res.m.a;
    //if refracted, ray has changed media
    if (ror < 0.) r.m = res.m;
}

void main() {
	//begin
    hit res;
    //final color
    vec3 final = vec3(0.);
    vec2 _uv, aspect = vec2(iResolution.x / iResolution.y, 1.);

	//scene
#ifdef SHOWSEG
	if (uv.y < 0.9) {
#endif
#ifdef MULTI
		for (int j = 0; j < MULTI; j++) {
			_uv = uv * aspect + (vec2(bell12(uv + float(j) - iGlobalTime*sml), bell12(uv - float(j) + iGlobalTime*sml)) * 2. - 1.) / iResolution.xy;
#else
			_uv = uv * aspect + (vec2(bell12(uv - iGlobalTime*sml), bell12(uv + iGlobalTime*sml)) * 2. - 1.) / iResolution.xy;
#endif

			mat startMat = _air;
			//start the ray at the camera
			ray r = ray(loc, normalize(rotateXY(normalize(vec3(_uv, 1.)), orient.xy)), vec3(EXPOSURE), startMat);
			//bounce around a few times
			for (int i = 0; i < BOUNCE_PXL; i++) {
				//new random seed  
#ifdef MULTI
				float rnd = float(iFrame*MULTI - i*BOUNCE_PXL + j)*sml;
#else
				float rnd = float(iFrame*BOUNCE_PXL - i)*sml;
#endif        
				//trace scene
				res = traceScene(r);
				//chance the medium refracts the light
				float rft = hash13(r.d + rnd);
				//probability increases with distance but tapers off
				if (rft < r.m.u.x * sqrt(res.d)) {
					//add emmisive color
	    			final += r.c * r.m.e * getEmission(r.d, res.n, r.m);
	    			//use brdf or btdf to scatter
					brdf(r, res, rnd);
					//burn a loop instead of doing another trace here for multi scattering
					continue;
				}
				//if collision
				if (res.d < zfar) {
					//update dynamic materials
					updateMaterials(res, r);
					//add colors
					final += r.c * res.m.e * getEmission(r.d, res.n, res.m);
					//modify ray properties
					bsdf(r, res, rnd);
				}
			}
#ifdef MULTI
		}
#endif
		final = pow(final, vec3(GAMMA));
#ifdef SHOWSEG
		//crosshair
        if (length(uv)<0.01) final = 1.-final;
#endif
		//continuous
		if (iMouse.z > 0. || length(vel) > sml) {
			fragColor = vec4(final, 1.);
			return;
		}
		vec4 lastFrame = texture2D(iChannel0, uv * .5 + .5);
		fragColor = vec4(final + lastFrame.rgb, lastFrame.a + 1.);
#ifdef SHOWSEG
	} else {
        //start the ray
        ray r = ray(loc, normalize(rotateXY(normalize(vec3(0., 0., 1.)), orient.xy)), vec3(1.), _air);
        //get segments of shapes
        seg a = tBox(r, box1), b = tBox(r, box2), c = tBox(r, box3);
        seg[] s = seg[6](nullSeg, nullSeg, nullSeg, nullSeg, nullSeg, nullSeg);
        tDiff(a,b, s[0], s[1]);
        tDiff(s[0],c, s[2], s[3]);
        tDiff(s[1],c, s[4], s[5]);
        //choose lowest non negitive point
        vec2 g = lt(lt(s[2], lt(s[3], lt(s[4], s[5]))));
        //choose lowest non negitive point
        float x = (uv.x+.25)*15., sx = .1*step(.95,fract(x)), uy = uv.y*2.-1.8;
		vec3 bgCol = (x < 0.)? vec3(.1+sx,0.,0.)+uy: vec3(.1)+sx+uy;
        vec3 fgcol = vec3(.1,vec2(step(0.,x)*.1));
        //render segments
        fragColor = vec4(bgCol + (
            inRange(s[0].t.x, x, s[0].t.y)+
            inRange(s[1].t.x, x, s[1].t.y)+
			inRange(s[2].t.x, x, s[2].t.y)+
            inRange(s[3].t.x, x, s[3].t.y)+
            inRange(s[4].t.x, x, s[4].t.y)+
            inRange(s[5].t.x, x, s[5].t.y)
        )*fgcol, 1.);
    }
#endif
}
