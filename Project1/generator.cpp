

#include <cmath>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>


constexpr double PI = 3.14159265358979323846;
inline constexpr double deg2rad(double deg) noexcept {
    return deg * PI / 180.0;
}


// A helper to generate uniformly distributed integers in [a, b].
static int randInt(int a, int b) {
    static std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

// A helper to generate uniformly distributed doubles in [a, b].
static double randDouble(double a, double b) {
    static std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

// A structure representing a question type and a function capable of
// constructing a question string given a difficulty level.
struct QuestionType {
    std::string name;
    std::function<std::string(int difficulty)> generator;
};

// A structure representing a category of questions, each with its own name
// and a list of question types.
struct Category {
    std::string name;
    std::vector<QuestionType> types;
};

// Forward declarations of generator functions.
static std::string generateRotationQuestion(int difficulty);
static std::string generateTranslationQuestion(int difficulty);
static std::string generateScalingQuestion(int difficulty);
static std::string generateColourQuestion(int difficulty);
static std::string generateShadingQuestion(int difficulty);
static std::string generateVectorOpsQuestion(int difficulty);
static std::string generateCompositionQuestion(int difficulty);
static std::string generateRayPlaneQuestion(int difficulty);
static std::string generateRaySphereQuestion(int difficulty);
static std::string generateMatrixFillQuestion(int difficulty);
static std::string generateProjectionQuestion(int difficulty);
static std::string generateColorConvertQuestion(int difficulty);


/*
 * Generate a rotation question.  The difficulty level controls the set
 * of permissible rotation angles and axes.  Difficulty 1 uses only 90° and
 * 180° rotations around the principal axes; difficulty 2 uses common angles
 * like 45° or 60°; difficulty 3 uses unusual integer angles (excluding those
 * in the easier categories); difficulty 4 may include an arbitrary axis
 * defined by a random unit vector and a non‑integer angle.
 */
static std::string generateRotationQuestion(int difficulty) {
    // Choose a base vector to rotate.
    std::vector<std::tuple<int, int, int>> vectors = {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}
    };
    auto [vx, vy, vz] = vectors[randInt(0, static_cast<int>(vectors.size()) - 1)];

    std::string axisDesc;
    double angle = 0.0;
    if (difficulty == 1) {
        angle = (randInt(0, 1) == 0 ? 90.0 : 180.0);
        int axis = randInt(0, 2);
        axisDesc = (axis == 0 ? "the X-axis" : (axis == 1 ? "the Y-axis" : "the Z-axis"));
    } else if (difficulty == 2) {
        angle = (randInt(0, 1) == 0 ? 45.0 : 60.0);
        int axis = randInt(0, 2);
        axisDesc = (axis == 0 ? "the X-axis" : (axis == 1 ? "the Y-axis" : "the Z-axis"));
    } else if (difficulty == 3) {
        int candidate;
        do {
            candidate = randInt(10, 80);
        } while (candidate % 15 == 0 || candidate == 45 || candidate == 60);
        angle = static_cast<double>(candidate);
        int axis = randInt(0, 2);
        axisDesc = (axis == 0 ? "the X-axis" : (axis == 1 ? "the Y-axis" : "the Z-axis"));
    } else {
        double ax = randDouble(-1.0, 1.0);
        double ay = randDouble(-1.0, 1.0);
        double az = randDouble(-1.0, 1.0);
        double length = std::sqrt(ax * ax + ay * ay + az * az);
        if (length == 0.0) length = 1.0;
        ax /= length;
        ay /= length;
        az /= length;
        angle = randDouble(5.0, 85.0);
        std::ostringstream oss;
        oss << "an arbitrary unit axis (" << std::fixed << std::setprecision(2)
            << ax << ", " << ay << ", " << az << ")";
        axisDesc = oss.str();
    }
    std::ostringstream q;
    q << "Given vector (" << vx << ", " << vy << ", " << vz
      << "), rotate it about " << axisDesc << " by " << angle
      << " degrees. Provide the coordinates of the rotated vector.";
    return q.str();
}

/*
 * Generate a translation question.
 */
static std::string generateTranslationQuestion(int difficulty) {
    std::ostringstream q;
    if (difficulty == 1) {
        int x = randInt(-3, 3);
        int y = randInt(-3, 3);
        int z = randInt(-3, 3);
        q << "Translate point (2, 1, 0) by vector (" << x << ", " << y << ", " << z
          << "). What are the new coordinates?";
    } else if (difficulty == 2) {
        double dx = randDouble(-5.0, 5.0);
        double dy = randDouble(-5.0, 5.0);
        double dz = randDouble(-5.0, 5.0);
        q << std::fixed << std::setprecision(2)
          << "Translate point (1.5, -2.0, 3.0) by vector (" << dx << ", " << dy << ", "
          << dz << "). What are the new coordinates?";
    } else if (difficulty == 3) {
        int x1 = randInt(-2, 2), y1 = randInt(-2, 2), z1 = randInt(-2, 2);
        int x2 = randInt(-2, 2), y2 = randInt(-2, 2), z2 = randInt(-2, 2);
        q << "Point (-1, 4, 0) is first translated by vector (" << x1 << ", " << y1 << ", " << z1
          << ") and then by vector (" << x2 << ", " << y2 << ", " << z2
          << "). What is its final position?";
    } else {
        int x_local = randInt(-3, 3), y_local = randInt(-3, 3), z_local = randInt(-3, 3);
        int frameRotation = randInt(1, 3) * 90;
        q << "In a local coordinate system, translate point (0, 2, 1) by vector (" << x_local << ", "
          << y_local << ", " << z_local << "). The local coordinate system is rotated "
          << frameRotation << " degrees about the Z-axis relative to the world coordinate system. Find the coordinates of the translated point in the world coordinate system.";
    }
    return q.str();
}

/*
 * Generate a scaling question.
 */
static std::string generateScalingQuestion(int difficulty) {
    std::ostringstream q;
    if (difficulty == 1) {
        int s = randInt(2, 5);
        q << "Scale point (2, -1, 3) uniformly by a factor of " << s << ". What are the scaled coordinates?";
    } else if (difficulty == 2) {
        int sx = randInt(2, 4);
        int sy = randInt(2, 4);
        int sz = randInt(2, 4);
        q << "Scale point (-1, 3, 2) along the X, Y and Z axes by factors (" << sx << ", "
          << sy << ", " << sz << "). What are the new coordinates?";
    } else if (difficulty == 3) {
        int sx = randInt(-3, -1);
        int sy = randInt(2, 5);
        int sz = randInt(2, 5);
        q << "Scale point (0, 2, -2) along the X, Y and Z axes by factors (" << sx << ", "
          << sy << ", " << sz << "). The negative factor along X indicates a mirror inversion. What is the result?";
    } else {
        int s1x = randInt(2, 4), s1y = randInt(2, 4), s1z = randInt(2, 4);
        int s2x = randInt(2, 4), s2y = randInt(2, 4), s2z = randInt(2, 4);
        q << "Point (-2, 1, 1) is first scaled by factors (" << s1x << ", " << s1y << ", " << s1z
          << ") and then by factors (" << s2x << ", " << s2y << ", " << s2z
          << "). What are the coordinates after the combined scaling?";
    }
    return q.str();
}

/*
 * Generate a basic colour theory question.
 */
static std::string generateColourQuestion(int difficulty) {
    std::ostringstream q;
    q.setf(std::ios::fixed);

    if (difficulty == 1) {
        int variant = randInt(0, 3);
        switch (variant) {
        case 0:
            q << "Briefly compare additive (RGB) vs subtractive (CMY/CMYK) colour models and give one typical use case for each.";
            break;
        case 1:
            q << "What are the additive primaries and secondary colours in the RGB model, and how are secondary colours formed?";
            break;
        case 2:
            q << "Explain gamma correction: why are display signals gamma-encoded, and why must values often be linearized before computation?";
            break;
        default:
            q << "Define luminance versus perceived brightness. Which one is physically defined and which one is perceptual?";
            break;
        }
    }
    else if (difficulty == 2) {
        int variant = randInt(0, 3);
        switch (variant) {
        case 0: {
            int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
            q << "Convert the RGB value (" << r << ", " << g << ", " << b
                << ") to HSV. Report H in degrees, S and V in [0,1].";
            break;
        }
        case 1: {
            int h = randInt(0, 359);
            double s = randDouble(0.2, 1.0);
            double v = randDouble(0.2, 1.0);
            q << std::setprecision(3)
                << "Convert HSV(" << h << "°, " << s << ", " << v
                << ") to RGB (0–255).";
            break;
        }
        case 2: {
            int r8 = randInt(0, 255);
            q << "Assume sRGB transfer. Given one channel value R8=" << r8
                << " in [0,255], normalize to [0,1], then linearize to get R_lin using the standard sRGB piecewise gamma. Report R_lin (4 decimals).";
            break;
        }
        default: {
            int r1 = randInt(0, 255), g1 = randInt(0, 255), b1 = randInt(0, 255);
            int r2 = randInt(0, 255), g2 = randInt(0, 255), b2 = randInt(0, 255);
            double a = randDouble(0.2, 0.8);
            q << std::setprecision(3)
                << "Linear-light alpha compositing (no premultiplied alpha): blend foreground C1=("
                << r1 << "," << g1 << "," << b1 << ") over background C2=("
                << r2 << "," << g2 << "," << b2 << ") with alpha=" << a
                << ". Steps: normalize, linearize (sRGB), compute C = a*C1 + (1-a)*C2, then re-encode to sRGB and scale to 0–255. Report blended RGB.";
            break;
        }
        }
    }
    else if (difficulty == 3) {
        int variant = randInt(0, 3);
        switch (variant) {
        case 0: {
            int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
            q << "Convert sRGB(" << r << "," << g << "," << b
                << ") to CIE XYZ under D65: normalize, linearize (sRGB), then multiply by the standard sRGB→XYZ matrix. Report X, Y, Z (3 decimals).";
            break;
        }
        case 1: {
            int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
            q << "Compute relative luminance Y for sRGB(" << r << "," << g << "," << b
                << "). Steps: normalize, linearize (sRGB), then Y = 0.2126 R + 0.7152 G + 0.0722 B. Report Y (3 decimals).";
            break;
        }
        case 2: {
            int X = randInt(10, 95), Y = randInt(10, 95), Z = randInt(10, 95);
            q << "Given XYZ(" << X << "," << Y << "," << Z
                << ") under D65/2°, convert to CIE L*a*b* using reference white (Xn=95.047, Yn=100.000, Zn=108.883). Report L*, a*, b* (1 decimal).";
            break;
        }
        default:
            q << "Explain chromatic adaptation (e.g., Bradford): what problem does it solve and at which step in a colour management pipeline it is applied?";
            break;
        }
    }
    else {
        int variant = randInt(0, 3);
        switch (variant) {
        case 0: {
            // ΔE*ab between two Lab colours
            int L1 = randInt(20, 90), a1 = randInt(-40, 40), b1 = randInt(-40, 40);
            int L2 = randInt(20, 90), a2 = randInt(-40, 40), b2 = randInt(-40, 40);
            q << "Given two Lab colours L1*=" << L1 << ", a1*=" << a1 << ", b1*=" << b1
                << " and L2*=" << L2 << ", a2*=" << a2 << ", b2*=" << b2
                << ", compute ΔE*ab = sqrt((ΔL*)^2 + (Δa*)^2 + (Δb*)^2). Report ΔE (2 decimals).";
            break;
        }
        case 1: {
            int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
            q << "Pipeline question: For sRGB(" << r << "," << g << "," << b
                << "), (1) linearize; (2) compute relative luminance Y; (3) decide which of two colours appears brighter by Y (you may compare against Y of sRGB(255,255,255)). Explain each step briefly.";
            break;
        }
        case 2: {
            q << "Why can’t a display with fixed RGB primaries cover the entire CIE xy chromaticity diagram? Discuss in terms of primary chromaticities, gamut triangle, and volume in a 3D colour space.";
            break;
        }
        default: {
            // HSV round trip: RGB -> HSV -> RGB to test precision and understanding
            int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
            q << "Round trip: Convert RGB(" << r << "," << g << "," << b
                << ") to HSV, then convert back to RGB and explain under what conditions small numeric differences may appear. Provide both intermediate HSV and final RGB (0–255).";
            break;
        }
        }
    }
    return q.str();
}

/*
 * Generate a shading question.
 */
static std::string generateShadingQuestion(int difficulty) {
    std::ostringstream q;
    if (difficulty == 1) {
        q << "In the Phong lighting model, what do the ambient, diffuse and specular components represent?";
    } else if (difficulty == 2) {
        q << "Given a surface normal (0, 0, 1), light direction (0, 1, 1), view direction (0, 0, 1),";
        q << " reflection coefficients ka=0.1, kd=0.7, ks=0.2, shininess exponent n=10 and white light intensity (1,1,1).";
        q << " Use the Phong model to compute the colour of the pixel (ignore distance attenuation).";
    } else if (difficulty == 3) {
        q << "Explain the difference between the Blinn–Phong model and the original Phong model, and why the half-vector improves computational efficiency.";
    } else {
        q << "To implement physically based rendering (PBR), which BRDF distribution functions are required? Explain the roles of the microfacet distribution and the geometry term in the Cook–Torrance model.";
    }
    return q.str();
}

static std::string generateVectorOpsQuestion(int difficulty) {
    std::ostringstream q;
    if (difficulty == 1) {
        int ax = randInt(-5, 5), ay = randInt(-5, 5), az = randInt(-5, 5);
        int bx = randInt(-5, 5), by = randInt(-5, 5), bz = randInt(-5, 5);
        if (randInt(0, 1) == 0) {
            q << "Given A(" << ax << "," << ay << "," << az << ") and B(" << bx << "," << by << "," << bz
                << "), compute A + B.";
        }
        else {
            q << "Given P(" << ax << "," << ay << "," << az << "), compute its distance to the origin.";
        }
    }
    else if (difficulty == 2) {
        int ax = randInt(-4, 4), ay = randInt(-4, 4), az = randInt(-4, 4);
        int bx = randInt(-4, 4), by = randInt(-4, 4), bz = randInt(-4, 4);
        if (ax == 0 && ay == 0 && az == 0) ax = 1; if (bx == 0 && by == 0 && bz == 0) bx = 1;
        if (randInt(0, 1) == 0) {
            q << "Given A(" << ax << "," << ay << "," << az << ") and B(" << bx << "," << by << "," << bz
                << "), compute the dot product A·B and the angle between them (degrees).";
        }
        else {
            q << "Given A(" << ax << "," << ay << "," << az << ") and B(" << bx << "," << by << "," << bz
                << "), determine whether they are orthogonal.";
        }
    }
    else if (difficulty == 3) {
        int ax = randInt(-3, 3), ay = randInt(-3, 3), az = randInt(-3, 3);
        int bx = randInt(-3, 3), by = randInt(-3, 3), bz = randInt(-3, 3);
        if (ax == 0 && ay == 0 && az == 0) ax = 1; if (bx == 0 && by == 0 && bz == 0) bx = 1;
        if (randInt(0, 1) == 0) {
            q << "Compute A×B for A(" << ax << "," << ay << "," << az << ") and B(" << bx << "," << by << "," << bz
                << "), and give the unit normal.";
        }
        else {
            int px = randInt(-3, 3), py = randInt(-3, 3), pz = randInt(-3, 3);
            q << "Project P(" << px << "," << py << "," << pz << ") onto the direction of A("
                << ax << "," << ay << "," << az << "). Provide the projected vector.";
        }
    }
    else {
        // 简单的正交基变换（给正交矩阵，算坐标）
        // 生成一个绕Z轴的正交旋转矩阵，避免数值坑
        int deg = (randInt(0, 3) + 1) * 30; // 30,60,90,120
        double c = std::cos(deg2rad(deg)), s = std::sin(deg2rad(deg));
        int vx = randInt(-3, 3), vy = randInt(-3, 3), vz = randInt(-3, 3);
        q << "Let R be a rotation about the Z-axis by " << deg << " degrees. "
            << "Given v(" << vx << "," << vy << "," << vz << "), compute R·v (coordinates in the rotated basis).";
    }
    return q.str();
}

static std::string generateCompositionQuestion(int difficulty) {
    std::ostringstream q;
    auto pickAngle = [&]() { int a[] = { 90,180,45,60 }; return a[randInt(0, 3)]; };
    if (difficulty == 1) {
        int tx = randInt(-3, 3), ty = randInt(-3, 3), tz = randInt(-3, 3);
        int ang = pickAngle(); char axis = "XYZ"[randInt(0, 2)];
        q << "Start from P(1,2,3). Apply translation T(" << tx << "," << ty << "," << tz
            << ") and then rotation R" << axis << "(" << ang << "°). Provide final coordinates.";
    }
    else if (difficulty == 2) {
        int tx = randInt(-3, 3), ty = randInt(-3, 3), tz = randInt(-3, 3);
        int ang = pickAngle(); char axis = "XYZ"[randInt(0, 2)];
        q << "For P(1,2,3), compare the results of R" << axis << "(" << ang << "°)∘T("
            << tx << "," << ty << "," << tz << ") vs T(" << tx << "," << ty << "," << tz << ")∘R" << axis
            << "(" << ang << "°). Explain whether they are identical and give both coordinates.";
    }
    else if (difficulty == 3) {
        int tx = randInt(-2, 2), ty = randInt(-2, 2), tz = randInt(-2, 2);
        int sx = randInt(-3, -1), sy = randInt(2, 4), sz = randInt(2, 4);
        int ang = pickAngle(); char axis = "XYZ"[randInt(0, 2)];
        q << "Given P(-1,0,2), apply S(" << sx << "," << sy << "," << sz << ") then R" << axis << "("
            << ang << "°) then T(" << tx << "," << ty << "," << tz << "). Provide final coordinates.";
    }
    else {
        int tx = randInt(-2, 2), ty = randInt(-2, 2), tz = randInt(-2, 2);
        int ang = pickAngle(); char axis = "XYZ"[randInt(0, 2)];
        q << "You are told that Q is obtained from P by T(" << tx << "," << ty << "," << tz
            << ") after R" << axis << "(" << ang << "°). Write the inverse sequence to recover P from Q "
            << "and compute P when Q(3,-1,4) is given.";
    }
    return q.str();
}

static std::string generateRayPlaneQuestion(int difficulty) {
    std::ostringstream q;
    // 射线 r(t) = O + tD
    int ox = randInt(-3, 3), oy = randInt(-3, 3), oz = randInt(0, 3);
    int dx = randInt(-2, 2), dy = randInt(-2, 2), dz = randInt(1, 3); 
    if (dx == 0 && dy == 0 && dz == 0) dz = 1;

    if (difficulty <= 2) {
        // 平面：n·x + d = 0  或  (x - P0)·n = 0
        int nx = 0, ny = 0, nz = 1; int d = -randInt(1, 4); // z = const
        if (difficulty == 2) {
            // 轻微变化法向
            int opts[3][3] = { {0,1,1},{1,0,1},{1,1,0} };
            int k = randInt(0, 2); nx = opts[k][0]; ny = opts[k][1]; nz = opts[k][2];
            if (nx == 0 && ny == 0 && nz == 0) nz = 1;
            d = -randInt(1, 4);
        }
        q << "Ray r(t)=O+tD with O(" << ox << "," << oy << "," << oz << ") and D("
            << dx << "," << dy << "," << dz << "). Plane: n·x + d = 0 with n("
            << nx << "," << ny << "," << nz << ") and d=" << d << ". Compute intersection t (if any) and point.";
    }
    else {
        // 点法式： (x - P0)·n = 0
        int p0x = randInt(-2, 2), p0y = randInt(-2, 2), p0z = randInt(1, 5);
        int nx = randInt(-1, 1), ny = randInt(-1, 1), nz = randInt(0, 1) + 1; // 让nz>=1
        q << "Ray r(t)=O+tD with O(" << ox << "," << oy << "," << oz << ") and D("
            << dx << "," << dy << "," << dz << "). Plane through P0(" << p0x << "," << p0y << "," << p0z
            << ") with normal n(" << nx << "," << ny << "," << nz << "). Compute t and intersection point; "
            << "also give the unit normal at the hit.";
    }
    return q.str();
}

static std::string generateRaySphereQuestion(int difficulty) {
    std::ostringstream q;
    int ox = randInt(-2, 2), oy = randInt(-2, 2), oz = randInt(-2, 2);
    int dx = randInt(-2, 2), dy = randInt(-2, 2), dz = randInt(1, 3); if (dx == 0 && dy == 0 && dz == 0) dz = 1;
    if (difficulty <= 2) {
        // 球心在原点
        int r = (difficulty == 1 ? 1 : randInt(1, 3));
        q << "Ray r(t)=O+tD with O(" << ox << "," << oy << "," << oz << ") and D("
            << dx << "," << dy << "," << dz << "). Sphere centered at origin with radius " << r
            << ". Compute intersection t value(s) and hit point(s).";
    }
    else {
        int cx = randInt(-2, 2), cy = randInt(-2, 2), cz = randInt(-2, 2), r = randInt(1, 3);
        if (cx == 0 && cy == 0 && cz == 0) cx = 1;
        q << "Ray r(t)=O+tD with O(" << ox << "," << oy << "," << oz << ") and D("
            << dx << "," << dy << "," << dz << "). Sphere center C(" << cx << "," << cy << "," << cz
            << ") with radius " << r << ". Compute number of intersections, t value(s), "
            << "and the surface normal at the first hit.";
    }
    return q.str();
}

static std::string generateMatrixFillQuestion(int difficulty) {
    std::ostringstream q;
    q.setf(std::ios::fixed); q << std::setprecision(4);

    if (difficulty == 1) {
     
        int kind = randInt(0, 1); // 0: translate, 1: scale
        if (kind == 0) {
            int tx = randInt(-3, 3), ty = randInt(-3, 3), tz = randInt(-3, 3);
            q << "Fill the 4x4 homogeneous translation matrix T(tx,ty,tz) for tx=" << tx
                << ", ty=" << ty << ", tz=" << tz << ". (Row-major or column-major is fine, "
                << "but be consistent and indicate where the translation terms go.)";
        }
        else {
            int sx = randInt(1, 4), sy = randInt(1, 4), sz = randInt(1, 4);
            q << "Fill the 4x4 scaling matrix S(sx,sy,sz) for sx=" << sx << ", sy=" << sy
                << ", sz=" << sz << ".";
        }
    }
    else if (difficulty == 2) {
       
        int sx = randInt(1, 4), sy = randInt(1, 4), sz = randInt(1, 4);
        int tx = randInt(-3, 3), ty = randInt(-3, 3), tz = randInt(-3, 3);
        int px = randInt(-2, 2), py = randInt(-2, 2), pz = randInt(-2, 2);
        int order = randInt(0, 1); // 0: T*S, 1: S*T
        q << "Let S=diag(" << sx << "," << sy << "," << sz << ",1) and T be translation by ("
            << tx << "," << ty << "," << tz << "). For point P(" << px << "," << py << "," << pz
            << ",1), compute (" << (order ? "S*T" : "T*S") << ")*P in homogeneous coordinates and give XYZ.";
    }
    else if (difficulty == 3) {
        
        int ang = (randInt(0, 3) + 1) * 30; // 30/60/90/120
        int px = randInt(-3, 3), py = randInt(-3, 3), pz = randInt(-1, 1);
        q << "Write the 4x4 rotation matrix Rz(" << ang << " deg) about Z, then apply it to P("
            << px << "," << py << "," << pz << ",1). Provide the rotated XYZ coordinates.";
    }
    else {
       
        int ang = (randInt(0, 3) + 1) * 45; // 45/90/135/180
        int tx = randInt(-3, 3), ty = randInt(-3, 3), tz = randInt(-3, 3);
        int qx = randInt(-3, 3), qy = randInt(-3, 3), qz = randInt(-3, 3);
        q << "A rigid transform M = T(" << tx << "," << ty << "," << tz << ") * Rz(" << ang
            << " deg). (Column-vector convention.) Write M^{-1} explicitly and compute P = M^{-1}*Q "
            << "for Q(" << qx << "," << qy << "," << qz << ",1). Give XYZ of P.";
    }
    return q.str();
}

static std::string generateProjectionQuestion(int difficulty) {
    std::ostringstream q;
    q.setf(std::ios::fixed); q << std::setprecision(4);

    if (difficulty == 1) {
        
        int fov = (randInt(0, 4) + 5) * 6; // 30,36,42,48,54
        double aspect = (randInt(0, 1) ? 16.0 / 9.0 : 4.0 / 3.0);
        double n = 0.1 * (randInt(1, 3));  // 0.1 / 0.2 / 0.3
        double f = 10.0 + randInt(0, 5) * 5; // 10,15,20,25,30,35
        q << "For a right-handed perspective projection with vertical FOV = " << fov << " deg, "
            << "aspect = " << aspect << ", near = " << n << ", far = " << f << ", write the 4x4 projection "
            << "matrix entries (m00, m11, m22, m23, m32). State the sign convention you use.";
    }
    else if (difficulty == 2) {
        
        int x = randInt(-3, 3), y = randInt(-3, 3); double z = -(randInt(1, 9)); // 负Z
        int fov = 45; double aspect = 16.0 / 9.0; double n = 0.1, f = 50.0;
        q << "Using a standard perspective matrix with FOV=45 deg, aspect=16/9, near=0.1, far=50, "
            << "project point P(" << x << "," << y << "," << z << ") to NDC (x_ndc,y_ndc,z_ndc). "
            << "Assume column-major, right-handed view looking down -Z.";
    }
    else if (difficulty == 3) {
       
        int l = -randInt(1, 4), r = randInt(2, 6), b = -randInt(1, 4), t = randInt(2, 6);
        double n = 0.1, f = 20.0;
        q << "Write the orthographic projection matrix for left=" << l << ", right=" << r
            << ", bottom=" << b << ", top=" << t << ", near=" << n << ", far=" << f
            << ". Then map P(" << randInt(l, r) << "," << randInt(b, t) << "," << -randInt(1, 10)
            << ") to NDC.";
    }
    else {
       
        int px = randInt(-3, 3), py = randInt(-3, 3), pz = -randInt(3, 10);
        int cx = randInt(-2, 2), cy = randInt(-2, 2), cz = randInt(1, 3); // camera position
        int yaw = (randInt(0, 7)) * 15;
        q << "Camera at C(" << cx << "," << cy << "," << cz << ") looking roughly along -Z with yaw " << yaw
            << " deg (no pitch/roll). Construct the view matrix V and perspective matrix P "
            << "(use FOV=60 deg, aspect=16/9, near=0.1, far=100). "
            << "Compute clip coords of point P(" << px << "," << py << "," << pz << ").";
    }
    return q.str();
}

static std::string generateColorConvertQuestion(int difficulty) {
    std::ostringstream q;
    q.setf(std::ios::fixed);

    if (difficulty == 1) {
     
        int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
        q << "Convert RGB(" << r << "," << g << "," << b << ") to HSV. Report H in degrees and S,V in [0,1]. "
            << "Assume sRGB primaries but ignore gamma (i.e., treat as already linear for HSV).";
    }
    else if (difficulty == 2) {
     
        int r = randInt(0, 255), g = randInt(0, 255), b = randInt(0, 255);
        q << "Convert sRGB(" << r << "," << g << "," << b << ") to CIE XYZ under D65. "
            << "Steps: normalize to [0,1], linearize sRGB (piecewise gamma), then multiply by the "
            << "standard sRGB-to-XYZ matrix. Give X,Y,Z (3 decimals).";
    }
    else if (difficulty == 3) {
       
        int X = randInt(10, 95), Y = randInt(10, 95), Z = randInt(10, 95);
        q << "Given XYZ(" << X << "," << Y << "," << Z << ") under D65/2°, convert to CIE L*a*b*. "
            << "Use the standard reference white (Xn=95.047, Yn=100.000, Zn=108.883). "
            << "Report L*, a*, b* (1 decimal).";
    }
    else {
        
        double L = randDouble(0.05, 0.95); // linear value
        double gamma = (randInt(0, 1) ? 2.2 : 2.4);
        q << "A linear-light value is L=" << std::setprecision(4) << L
            << ". Apply display gamma " << gamma << " to obtain the encoded sRGB channel value E "
            << "(i.e., E=L^(1/gamma), ignore sRGB's piecewise segment). Then invert it back to "
            << "linear. Report both E and the recovered linear value.";
    }
    return q.str();
}


/*
 * Build the catalogue of categories and question types.
 */
static std::vector<Category> buildCatalogue() {
    std::vector<Category> catalogue;

    Category geom;
    // Use English names for categories and question types
    geom.name = "Geometric Transformations";
    geom.types.push_back({"Rotation", generateRotationQuestion});
    geom.types.push_back({"Translation", generateTranslationQuestion});
    geom.types.push_back({"Scaling", generateScalingQuestion});
    geom.types.push_back({ "Vector Operations", generateVectorOpsQuestion });
    geom.types.push_back({ "Composition Order", generateCompositionQuestion });
    geom.types.push_back({ "Matrix Fill / Apply", generateMatrixFillQuestion });
    catalogue.push_back(geom);

    Category colour;
    colour.name = "Color Theory";
    colour.types.push_back({"Color Models & Conversion", generateColourQuestion});
    colour.types.push_back({ "Color Conversion (Numeric)", generateColorConvertQuestion });
    catalogue.push_back(colour);

    Category shading;
    shading.name = "Lighting & Shading";
    shading.types.push_back({"Shading Models", generateShadingQuestion});
    catalogue.push_back(shading);

    Category ray;
    ray.name = "Ray Tracing (Core)";
    ray.types.push_back({ "Ray–Plane Intersection", generateRayPlaneQuestion });
    ray.types.push_back({ "Ray–Sphere Intersection", generateRaySphereQuestion });
    catalogue.push_back(ray);

    Category proj;
    proj.name = "Projection & Viewing";
    proj.types.push_back({ "Perspective/Ortho Projection", generateProjectionQuestion });
    catalogue.push_back(proj);


    return catalogue;
}

/*
 * Main entry point.  Presents an interactive menu for generating questions.
 */
int main() {
    std::vector<Category> catalogue = buildCatalogue();
    std::cout << "Welcome to the exam question generator.\n";
    bool running = true;
    while (running) {
        std::cout << "\nPlease choose a category:\n";
        for (std::size_t i = 0; i < catalogue.size(); ++i) {
            std::cout << "  " << (i + 1) << ": " << catalogue[i].name << "\n";
        }
        std::cout << "  0: Exit\n";
        std::cout << "Enter the number: ";
        int catIndex;
        std::cin >> catIndex;
        if (!std::cin.good()) {
            std::cin.clear(); std::cin.ignore(1024, '\n');
            std::cerr << "Invalid input, please try again.\n";
            continue;
        }
        if (catIndex == 0) {
            running = false;
            break;
        }
        if (catIndex < 1 || static_cast<std::size_t>(catIndex) > catalogue.size()) {
            std::cerr << "Invalid category number, please try again.\n";
            continue;
        }
        Category &cat = catalogue[catIndex - 1];
        std::cout << "\nIn category [" << cat.name << "] please choose a question type:\n";
        for (std::size_t i = 0; i < cat.types.size(); ++i) {
            std::cout << "  " << (i + 1) << ": " << cat.types[i].name << "\n";
        }
        std::cout << "Enter the number: ";
        int typeIndex;
        std::cin >> typeIndex;
        if (!std::cin.good()) {
            std::cin.clear(); std::cin.ignore(1024, '\n');
            std::cerr << "Invalid input, please try again.\n";
            continue;
        }
        if (typeIndex < 1 || static_cast<std::size_t>(typeIndex) > cat.types.size()) {
            std::cerr << "Invalid question type number, please try again.\n";
            continue;
        }
        QuestionType &qtype = cat.types[typeIndex - 1];
        std::cout << "\nChoose a difficulty level (1-4): ";
        int difficulty;
        std::cin >> difficulty;
        if (!std::cin.good() || difficulty < 1 || difficulty > 4) {
            std::cin.clear(); std::cin.ignore(1024, '\n');
            std::cerr << "Invalid difficulty level. Please enter 1-4.\n";
            continue;
        }
        std::cout << "Enter the number of questions to generate: ";
        int count;
        std::cin >> count;
        if (!std::cin.good() || count < 1) {
            std::cin.clear(); std::cin.ignore(1024, '\n');
            std::cerr << "The number of questions must be a positive integer.\n";
            continue;
        }
        std::cout << "\nGenerated questions:\n";
        for (int i = 0; i < count; ++i) {
            std::string question = qtype.generator(difficulty);
            std::cout << (i + 1) << ". " << question << "\n\n";
        }
    }
    std::cout << "Thank you for using the generator. Goodbye!\n";
    return 0;
}