#include "trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;
//#include <getopt.h>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
//#include <values.h>
#define MAX DBL_MAX
#endif

SlVector3 normalizeVector(SlVector3 a) {
    SlVector3 newVector;
    double num = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    newVector[0] = a[0] / num;
    newVector[1] = a[1] / num;
    newVector[2] = a[2] / num;
    return newVector;
}
// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3 &a, const SlVector3 &b, const SlVector3 &c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
        b[0] * (c[1] * a[2] - a[1] * c[2]) +
            c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) {return x*x;} 

bool Triangle::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {
    SlVector3 start = r.e;//starting point
    SlVector3 direction = r.d;//direction of ray
    //v1 v2 v3 is a,b,c
    // Step 1 Ray-triangle test
    SlVector3 firstLine, secondLine, thirdLine;
    SlVector3 betas = b - a;
    SlVector3 gammas = c - a;
    SlVector3 normal = normalizeVector(cross((b - a), (c - a)));
    double fConstant, sConstant, tConstant;//constants for the first, second and third lines;
    firstLine[0] = betas[0]; firstLine[1] = gammas[0]; firstLine[2] = 0 - direction[0];
    secondLine[0] = betas[1]; secondLine[1] = gammas[1]; secondLine[2] = 0 - direction[1];
    thirdLine[0] = betas[2]; thirdLine[1] = gammas[2]; thirdLine[2] = 0 - direction[2];
    fConstant = start[0] - a[0];
    sConstant = start[1] - a[1];
    tConstant = start[2] - a[2];
    //cout << firstLine[0] <<" " << firstLine[1] << " " << firstLine[2] << " " << fConstant << endl;
    //cout << secondLine[0] << " " << secondLine[1] << " " << secondLine[2] << " " << sConstant << endl;
    //cout << thirdLine[0] << " " << thirdLine[1] << " " << thirdLine[2] << " " << tConstant << endl;
    double xOne = firstLine[0] * secondLine[2] - secondLine[0] * firstLine[2];
    double yOne = firstLine[1] * secondLine[2] - secondLine[1] * firstLine[2];
    double constOne = secondLine[2] * fConstant - firstLine[2] * sConstant;
    double xTwo = firstLine[0] * thirdLine[2] - thirdLine[0] * firstLine[2];
    double yTwo = firstLine[1] * thirdLine[2] - thirdLine[1] * firstLine[2];
    double constTwo = thirdLine[2] * fConstant - firstLine[2] * tConstant;
    double beta = (constOne * yTwo - constTwo * yOne) / (xOne * yTwo - xTwo * yOne);
    double gamma = (constOne - xOne * beta) / yOne;
    double t = (fConstant - firstLine[0] * beta - firstLine[1] * gamma) / firstLine[2];
    //cout << beta << " " << gamma << " " << t << endl; //bgt
    if (beta + gamma <= 1 && t > t0 && t < t1) {
        if (hr.t != 0 && t < hr.t) {
            SlVector3 point = start + t * direction;
            hr.p = point;
            hr.t = t;
            hr.n = normal;
            return true;
        }
    }

    return false;
}

bool TrianglePatch::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    bool temp = Triangle::intersect(r,t0,t1,hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
        return true;
    }
    return false;
}


bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {

    // Step 1 Sphere-ray test
    SlVector3 startPoint = r.e;
    SlVector3 Direction = r.d;
    double constTerm = dot(startPoint,startPoint)+dot(c,c)-rad*rad-2*dot(startPoint,c);//c
    double oneTerm = 2 * (dot(startPoint, Direction) - dot(c, Direction));//b
    double secondTerm = dot(Direction, Direction);//a
    //cout << constTerm << " " << oneTerm << " " << secondTerm << " "<<endl;
    double delta = oneTerm * oneTerm - 4 * constTerm * secondTerm;
    if (delta > 0) {
        cout << delta << "delta" << endl;
    }
    if (delta < 0) {
        return false;
    }
    if (delta == 0) {
        double answer = (0 - oneTerm) / (2 * secondTerm);
        //cout << answer<<"=0ans"<<endl;
        if (answer > t0 && answer < t1) {
            if (hr.t != 0 && answer < hr.t) {
                SlVector3 point = startPoint + answer * Direction;
                hr.p = point;
                hr.t = answer;
                hr.n = normalizeVector(point - c);
                return true;
            }
        }
        else {
            return false;
        }

    }
    else {
        double answer1 = (0-sqrt(delta) - oneTerm) / (2 * secondTerm);
        double answer2 = (sqrt(delta) - oneTerm) / (2 * secondTerm);
        double answer;
        if (answer1 > t0 && answer1 < t1) {
            answer = answer1;
            //cout << answer << "ans1" << endl;
        }
        else if (answer2 > t0 && answer2 < t1) {
            answer = answer2;
            //cout << answer << "ans2" << endl;
        }
        
        else {
            return false;
        }
        if (hr.t != 0 && answer < hr.t) {
            SlVector3 point = startPoint + answer * Direction;
            hr.p = point;
            hr.t = answer;
            hr.n = normalizeVector(point - c);
            return true;
        }
    }
    return false;
}




Tracer::Tracer(const std::string &fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        getline(in, line);
        switch (line[0]) {
            case 'b': {
                std::stringstream ss(line);
                ss>>ch>>bcolor[0]>>bcolor[1]>>bcolor[2];
                break;
            }

            case 'v': {
                getline(in, line);
                std::string junk;
                std::stringstream fromss(line);
                fromss>>junk>>eye[0]>>eye[1]>>eye[2];

                getline(in, line);
                std::stringstream atss(line);
                atss>>junk>>at[0]>>at[1]>>at[2];

                getline(in, line);
                std::stringstream upss(line);
                upss>>junk>>up[0]>>up[1]>>up[2];

                getline(in, line);
                std::stringstream angless(line);
                angless>>junk>>angle;

                getline(in, line);
                std::stringstream hitherss(line);
                hitherss>>junk>>hither;

                getline(in, line);
                std::stringstream resolutionss(line);
                resolutionss>>junk>>res[0]>>res[1];
                break;
            }

            case 'p': {
                bool patch = false;
                std::stringstream ssn(line);
                unsigned int nverts;
                if (line[1] == 'p') {
                    patch = true;
                    ssn>>ch;
                }
                ssn>>ch>>nverts;
                std::vector<SlVector3> vertices;
                std::vector<SlVector3> normals;
                for (unsigned int i=0; i<nverts; i++) {
                    getline(in, line);
                    std::stringstream ss(line);
                    SlVector3 v,n;
                    if (patch) ss>>v[0]>>v[1]>>v[2]>>n[0]>>n[1]>>n[2];
                    else ss>>v[0]>>v[1]>>v[2];
                    vertices.push_back(v);
                    normals.push_back(n);
                }
                bool makeTriangles = false;
                if (vertices.size() == 3) {
                    if (patch) {
                        surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
                        normals [0], normals [1], normals [2]), fill));
                    } else {
                        surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                    }
                }
 else if (vertices.size() == 4) {
 SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
 SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
 SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
 SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
 if (dot(n0, n1) > 0 && dot(n0, n2) > 0 && dot(n0, n3) > 0) {
     makeTriangles = true;
     if (patch) {
         surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
             normals[0], normals[1], normals[2]), fill));
         surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3],
             normals[0], normals[2], normals[3]), fill));
     }
     else {
         surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
         surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
     }
 }
 if (!makeTriangles) {
     std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
 }
                }
                break;
            }

            case 's': {
                std::stringstream ss(line);
                SlVector3 c;
                double r;
                ss >> ch >> c[0] >> c[1] >> c[2] >> r;
                surfaces.push_back(std::pair<Surface*, Fill>(new Sphere(c, r), fill));
                break;
            }

            case 'f': {
                std::stringstream ss(line);
                ss >> ch >> fill.color[0] >> fill.color[1] >> fill.color[2] >> fill.kd >> fill.ks >> fill.shine >> fill.t >> fill.ior;
                break;
            }

            case 'l': {
                std::stringstream ss(line);
                Light l;
                ss >> ch >> l.p[0] >> l.p[1] >> l.p[2];
                if (!ss.eof()) {
                    ss >> l.c[0] >> l.c[1] >> l.c[2];
                    coloredlights = true;
                }
                lights.push_back(l);
                break;
            }

            default:
                break;
        }
    }
    if (!coloredlights) for (unsigned int i = 0; i < lights.size(); i++) lights[i].c = 1.0 / sqrt(lights.size());
    im = new SlVector3[res[0] * res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}

Tracer::~Tracer() {
    if (im) delete[] im;
    for (unsigned int i = 0; i < surfaces.size(); i++) delete surfaces[i].first;
}

SlVector3 Tracer::shade(const HitRecord& hr) const {
    if (color) return hr.f.color;

    SlVector3 color(0.0);
    HitRecord dummy;

    for (unsigned int i = 0; i < lights.size(); i++) {
        const Light& light = lights[i];
        bool shadow = false;

        // Step 3 Check for shadows here

        if (!shadow) {

            // Step 2 do shading here

        }
    }


    // Step 4 Add code for computing reflection color here

    // Step 5 Add code for computing refraction color here

    return color;
}
/*
SlVector3 Tracer::shade(const HitRecord& hr) const {
    if (color) return hr.f.color;

    SlVector3 color(0.0);
    HitRecord dummy;
    /*
    for (unsigned int i = 0; i < lights.size(); i++) {
        const Light& light = lights[i];
        bool shadow = false;

        // Step 3 Check for shadows here
        SlVector3 surfacetoLight = hr.p - light.p;
        Ray surfacetoLightRay = Ray(hr.p, normalizeVector(surfacetoLight),0);
        //to do: check intersection
        bool ifIntersect;
        if (ifIntersect) {
            shadow = true;
            return color;
        }
        if (!shadow) {

            // Step 2 do shading here

        }
    }


    // Step 4 Add code for computing reflection color here

    // Step 5 Add code for computing refraction color here
    
    return color;
}
*/
SlVector3 Tracer::trace(const Ray& r, double t0, double t1) const {
    HitRecord hr;
    SlVector3 color(bcolor);
    //std::vector<std::pair<Surface *, Fill> > surfaces;
    bool hit = false;
    for (int i = 0; i < surfaces.size(); i++) {
        bool currentHit = surfaces[i].first->intersect(r, t0, t1, hr);
        if (currentHit) {
            hit = true;
            cout << "yes";
        }
        // Step 1 See what a ray hits      
        
    }
    if (hit) {
        hr.raydepth += 1;
        //color = shade(hr);
        color = hr.f.color;
    }
    return color;
}

void Tracer::traceImage() {
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up,w);
    normalize(u);
    SlVector3 v = cross(w,u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((3.14/180.0) * (angle/2.0)) * d; //Originally M_Pi
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3 *pixel = im;

    for (unsigned int j=0; j<res[1]; j++) {
        for (unsigned int i=0; i<res[0]; i++, pixel++) {

            SlVector3 result(0.0,0.0,0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r-l)*(i+rx)/res[0];
                double y = b + (t-b)*(j+ry)/res[1];
                SlVector3 dir = -d * w + x * u + y * v;
	
                Ray r(eye, dir);
                normalize(r.d);

                result += trace(r, hither, MAX);

            }
            (*pixel) = result / samples;
        }
    }
}

void Tracer::writeImage(const std::string &fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out<<"P6"<<"\n"<<res[0]<<" "<<res[1]<<"\n"<<255<<"\n";
    SlVector3 *pixel = im;
    char val;
    for (unsigned int i=0; i<res[0]*res[1]; i++, pixel++) {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write (&val, sizeof(unsigned char));
    }
    out.close();
}


int main(int argc, char *argv[]) {
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    /*
    while ((c = getopt(argc, argv, "a:s:d:c")) != -1) {
        switch(c) {
            case 'a':
            aperture = atof(optarg);
            break;
            case 's':
            samples = atoi(optarg);
            break;
            case 'c':
            color = true;
            break;
            case 'd':
            maxraydepth = atoi(optarg);
            break;
            default:
            abort();
        }
    }
    
    
    if (argc-optind != 2) {
        std::cout<<"usage: trace [opts] input.nff output.ppm"<<std::endl;
        for (unsigned int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
        exit(0);
    }	
    
    
    SlVector3 center;
    center[0] = 3;
    center[1] = 3;
    center[2] = 0;
    SlVector3 e;
    SlVector3 d;
    e[0] = 0; e[1] = 0; e[2] = 0;
    d[0] = -1; d[1] = -1; d[2] = 0;
    Sphere testS(center,2);
    Ray testRay(e, d, 0);
    HitRecord hr;
    cout << testS.intersect(testRay, 0, 1000000, hr);
    
    SlVector3 point1; point1[0] = 3; point1[1] = 1; point1[2] = 1;
    SlVector3 point2; point2[0] = 7; point2[1] = 2; point2[2] = 2;
    SlVector3 point3; point3[0] = 4; point3[1] = 4; point3[2] = 3;
    SlVector3 e;
    SlVector3 d;
    e[0] = 0; e[1] = 1; e[2] = 0;
    d[0] = 1; d[1] = 0; d[2] = 0;

    Ray testRay(e, d, 0);
    HitRecord hr;
    Triangle testT(point1, point2, point3);
    cout << testT.intersect(testRay, 0, 10000, hr);
    */
    Tracer tracer(argv[1]);
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    tracer.traceImage();
    tracer.writeImage(argv[2]);
    
};
