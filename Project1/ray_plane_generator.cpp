//// ray_plane_generator.cpp -- MSVC-friendly, with D3 template + robust math + seed + export
//// Project: C/C++ -> Language -> /std:c++17
//
//#include <iostream>
//#include <iomanip>
//#include <sstream>
//#include <string>
//#include <vector>
//#include <array>
//#include <random>
//#include <ctime>
//#include <cmath>
//#include <utility>
//#include <functional>
//#include <limits>
//#include <fstream>   // export
//#include <optional>  // if you need optional later
//
//using namespace std;
//
//// ---------------- RNG & seed ----------------
//static mt19937 rng(static_cast<unsigned>(std::time(nullptr)));
//static void reseed(unsigned s) { rng.seed(s); }
//static int randi(int a, int b) { uniform_int_distribution<int> dist(a, b); return dist(rng); }
//static double randf(double a, double b) { uniform_real_distribution<double> dist(a, b); return dist(rng); }
//
//// ---------------- Config ----------------
//static constexpr double EPS = 1e-6;   // 统一容差
//
//// ---------------- Rubric & Difficulty ----------------
//struct Rubric {
//    int form = 0, generality = 0, branching = 0, steps = 0, downstream = 0;
//    int total() const { return form + generality + branching + steps + downstream; }
//};
//enum class Difficulty { D1 = 1, D2, D3, D4 };
//static pair<int, int> bandFor(Difficulty d) {
//    switch (d) {
//    case Difficulty::D1: return { 0,3 };
//    case Difficulty::D2: return { 4,6 };
//    case Difficulty::D3: return { 7,9 };
//    default:             return { 10,12 };
//    }
//}
//enum class VizLevel { None = 0, Qualitative = 1 }; // 预留
//
//// ---------------- Payload / Result ----------------
//struct RayPlaneParams {
//    array<double, 3> o{ 0,0,0 }, d{ 0,0,1 }; // ray
//    // plane (we will set one of the forms per template)
//    double zc = std::numeric_limits<double>::quiet_NaN(); // z=c
//    array<double, 3> n{ 0,0,1 };  // normal
//    double d0 = 0;             // n·x + d0 = 0
//    array<double, 3> p0{ 0,0,0 }; // point on plane (point-normal form)
//};
//
//struct Generated {
//    string prompt;
//    string answer;
//    string metadata;
//};
//
//// ---------------- small vector helpers ----------------
//static string vec3s(const array<double, 3>& v, int prec = 3) {
//    ostringstream os; os.setf(std::ios::fixed); os << setprecision(prec);
//    os << "(" << v[0] << "," << v[1] << "," << v[2] << ")"; return os.str();
//}
//static array<double, 3> add(const array<double, 3>& a, const array<double, 3>& b) {
//    return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
//}
//static array<double, 3> sub(const array<double, 3>& a, const array<double, 3>& b) {
//    return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
//}
//static array<double, 3> mul(const array<double, 3>& a, double t) {
//    return { a[0] * t, a[1] * t, a[2] * t };
//}
//static double dot(const array<double, 3>& a, const array<double, 3>& b) {
//    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//static array<double, 3> normalize(const array<double, 3>& v) {
//    double L = std::sqrt(std::max(1e-18, dot(v, v)));
//    return { v[0] / L, v[1] / L, v[2] / L };
//}
//
//// ---------------- robust math (统一求交) ----------------
//enum class HitState { Hit, Parallel, Coplanar, Behind, NoHit };
//
//struct HitResult {
//    HitState state{ HitState::NoHit };
//    double t = std::numeric_limits<double>::quiet_NaN();
//    array<double, 3> hit{ 0,0,0 };
//    array<double, 3> unitN{ 0,0,1 };
//};
//
//static HitResult intersect_general(const array<double, 3>& o, const array<double, 3>& d,
//    const array<double, 3>& n, double d0, double eps = EPS) {
//    HitResult R; R.unitN = normalize(n);
//    double denom = dot(n, d);
//    if (std::fabs(denom) < eps) {
//        double val = dot(n, o) + d0;
//        if (std::fabs(val) < eps) { R.state = HitState::Coplanar; }
//        else { R.state = HitState::Parallel; }
//        return R;
//    }
//    double numer = -(dot(n, o) + d0);
//    double t = numer / denom;
//    if (t < 0) { R.state = HitState::Behind; R.t = t; return R; }
//    R.state = HitState::Hit; R.t = t; R.hit = add(o, mul(d, t));
//    return R;
//}
//
//static HitResult intersect_pointnormal(const array<double, 3>& o, const array<double, 3>& d,
//    const array<double, 3>& p0, const array<double, 3>& n,
//    double eps = EPS) {
//    HitResult R; R.unitN = normalize(n);
//    double denom = dot(d, n);
//    if (std::fabs(denom) < eps) {
//        double rel = dot(sub(o, p0), n);
//        if (std::fabs(rel) < eps) { R.state = HitState::Coplanar; }
//        else { R.state = HitState::Parallel; }
//        return R;
//    }
//    double numer = dot(sub(p0, o), n);
//    double t = numer / denom;
//    if (t < 0) { R.state = HitState::Behind; R.t = t; return R; }
//    R.state = HitState::Hit; R.t = t; R.hit = add(o, mul(d, t));
//    return R;
//}
//
//// ---------------- parameter budgeter（避免近平行） ----------------
//static array<double, 3> sample_dir_not_parallel_to(const array<double, 3>& n, double eps = 0.1) {
//    // 随机直到 |dot(d,n)| >= eps
//    while (true) {
//        array<double, 3> d{ randf(-1,1), randf(-1,1), randf(-1,1) };
//        double L = std::sqrt(std::max(1e-12, dot(d, d)));
//        d = { d[0] / L, d[1] / L, d[2] / L };
//        if (std::fabs(dot(d, n)) >= eps) return d;
//    }
//}
//
//// ---------------- metadata ----------------
//static string mkMetaJSON(const string& name, const Rubric& r, Difficulty D,
//    const RayPlaneParams& P, VizLevel viz = VizLevel::None) {
//    ostringstream os; os.setf(std::ios::fixed); os << setprecision(6);
//    os << "{";
//    os << "\"type\":\"Ray-Plane\"";
//    os << ",\"template\":\"" << name << "\"";
//    os << ",\"rubric\":{\"form\":" << r.form << ",\"generality\":" << r.generality
//        << ",\"branching\":" << r.branching << ",\"steps\":" << r.steps
//        << ",\"downstream\":" << r.downstream << ",\"total\":" << r.total() << "}";
//    os << ",\"difficulty\":\"D" << int(D) << "\"";
//    os << ",\"epsilon\":" << EPS;
//    os << ",\"params\":{";
//    os << "\"o\":[" << P.o[0] << "," << P.o[1] << "," << P.o[2] << "],";
//    os << "\"d\":[" << P.d[0] << "," << P.d[1] << "," << P.d[2] << "],";
//    os << "\"zc\":" << (std::isnan(P.zc) ? -9999.0 : P.zc) << ",";
//    os << "\"n\":[" << P.n[0] << "," << P.n[1] << "," << P.n[2] << "],";
//    os << "\"d0\":" << P.d0 << ",";
//    os << "\"p0\":[" << P.p0[0] << "," << P.p0[1] << "," << P.p0[2] << "]}";
//    os << ",\"viz\":\"" << (viz == VizLevel::None ? "None" : "Qualitative") << "\"";
//    os << "}";
//    return os.str();
//}
//
//// ---------------- Template interface ----------------
//struct RayPlaneTemplate {
//    string name;
//    Rubric rubric;
//    function<Generated()> instantiate;
//};
//
//// D1: z=c, guaranteed hit; 要学生写出 t 公式与意义（理解优先）
//static RayPlaneTemplate T_D1_axisZ() {
//    RayPlaneTemplate T;
//    T.name = "AxisZ_known_hit_t_and_point";
//    T.rubric = { 0,0,0,0,0 }; // total 0 → D1
//    T.instantiate = [=]() {
//        RayPlaneParams P;
//        P.o = { double(randi(-2,2)), double(randi(-2,2)), double(randi(0,2)) };
//        P.d = { double(randi(-1,1)), double(randi(-1,1)), double(randi(1,3)) };
//        P.zc = double(randi(int(P.o[2]) + 1, int(P.o[2]) + 4));
//
//        // 等价转成一般式 n·x + d0 = 0
//        P.n = { 0,0,1 }; P.d0 = -P.zc;
//
//        auto R = intersect_general(P.o, P.d, P.n, P.d0, EPS);
//
//        ostringstream q, a; q.setf(std::ios::fixed); a.setf(std::ios::fixed);
//        q << "Ray r(t)=o+t d with o=" << vec3s(P.o) << ", d=" << vec3s(P.d)
//            << "; Plane z=" << P.zc << ".\n"
//            << "Tasks: (1) Write the formula for t when the plane is z=c, and explain why it follows from o_z + t d_z = c.\n"
//            << "(2) Using ε=" << EPS << ", compute t and the hit point.";
//
//        a << setprecision(6);
//        if (R.state == HitState::Hit) a << "t=" << R.t << ", hit=" << vec3s(R.hit, 3);
//        else a << "Unexpected state (should be hit in D1).";
//        return Generated{ q.str(), a.str(), mkMetaJSON(T.name, T.rubric, Difficulty::D1, P) };
//        };
//    return T;
//}
//
//// D2: n·x + d0 = 0（小整数），要求并解释“分母检查”和最近非负 t
//static RayPlaneTemplate T_D2_general_ndotx_plus_d() {
//    RayPlaneTemplate T;
//    T.name = "GeneralForm_smallInt_check_denominator";
//    T.rubric = { 1,1,1,1,1 }; // total 5 → D2
//    T.instantiate = [=]() {
//        RayPlaneParams P;
//        P.o = { double(randi(-2,2)), double(randi(-2,2)), double(randi(-2,2)) };
//        P.n = { double(randi(0,1)), double(randi(0,1)), 1.0 }; // 友好法向
//        P.d = sample_dir_not_parallel_to(P.n, 0.2);           // 避免近平行
//        P.d0 = -double(randi(1, 4));
//
//        auto R = intersect_general(P.o, P.d, P.n, P.d0, EPS);
//
//        ostringstream q, a; q.setf(std::ios::fixed); a.setf(std::ios::fixed);
//        q << "Ray r(t)=o+t d, o=" << vec3s(P.o) << ", d=" << vec3s(P.d)
//            << "; Plane n·x + d = 0 with n=" << vec3s(P.n) << ", d=" << P.d0 << ".\n"
//            << "Tasks: (1) Check parallelism by n·d≈0 with ε=" << EPS << "; explain why this matters.\n"
//            << "(2) If not parallel, compute the **nearest non-negative** t and the hit point.";
//
//        a << setprecision(6);
//        if (R.state == HitState::Parallel) a << "Parallel (n·d≈0). No unique intersection.";
//        else if (R.state == HitState::Behind) a << "t=" << R.t << " < 0 (behind origin). No forward hit.";
//        else if (R.state == HitState::Hit) a << "t=" << R.t << ", hit=" << vec3s(R.hit, 3);
//        else if (R.state == HitState::Coplanar) a << "Coplanar: ray lies in the plane.";
//        else a << "No hit.";
//        return Generated{ q.str(), a.str(), mkMetaJSON(T.name, T.rubric, Difficulty::D2, P) };
//        };
//    return T;
//}
//
//// D3: (x-p0)·n=0，需分类 + 最近非负 t（不要求单位法线）——偏“理解+判别”
//static RayPlaneTemplate T_D3_pointNormal_classify_nearestT() {
//    RayPlaneTemplate T;
//    T.name = "PointNormal_classify_nearestT_noUnitNormal";
//    T.rubric = { 2,2,2,1,1 }; // total 8 → D3
//    T.instantiate = [=]() {
//        RayPlaneParams P;
//        P.p0 = { double(randi(-2,2)), double(randi(-2,2)), double(randi(0,4)) };
//        P.n = { double(randi(0,1)), double(randi(1,2)), double(randi(1,2)) };
//        // 光线方向避免近平行，让学生能专注“分类与解释”
//        P.d = sample_dir_not_parallel_to(P.n, 0.15);
//        P.o = { double(randi(-2,2)), double(randi(-2,2)), double(randi(-2,2)) };
//
//        auto R = intersect_pointnormal(P.o, P.d, P.p0, P.n, EPS);
//
//        ostringstream q, a; q.setf(std::ios::fixed); a.setf(std::ios::fixed);
//        q << "Ray r(t)=o+t d, o=" << vec3s(P.o) << ", d=" << vec3s(P.d)
//            << "; Plane (x - p0)·n = 0 with p0=" << vec3s(P.p0) << ", n=" << vec3s(P.n) << ".\n"
//            << "Tasks: (1) Classify the case: parallel / coplanar / hit / no-forward-hit, using ε=" << EPS << ".\n"
//            << "(2) If hit, compute the **nearest non-negative** t and the intersection point.\n"
//            << "(3) Briefly explain the sign meaning of t (why t<0 is \"behind origin\").";
//
//        a << setprecision(6);
//        if (R.state == HitState::Parallel) a << "Parallel (|d·n|<ε).";
//        else if (R.state == HitState::Coplanar) a << "Coplanar: infinite intersections.";
//        else if (R.state == HitState::Behind) a << "t=" << R.t << " < 0 (behind origin).";
//        else if (R.state == HitState::Hit) a << "t=" << R.t << ", hit=" << vec3s(R.hit, 3);
//        else a << "No hit.";
//        return Generated{ q.str(), a.str(), mkMetaJSON(T.name, T.rubric, Difficulty::D3, P) };
//        };
//    return T;
//}
//
//// D4: (x-p0)·n=0，分类 + 最近非负 t + **单位法线**（整合更多要素）
//static RayPlaneTemplate T_D4_pointNormal_with_unitNormal() {
//    RayPlaneTemplate T;
//    T.name = "PointNormal_classify_and_unit_normal";
//    T.rubric = { 2,2,2,2,2 }; // total 10 → D4
//    T.instantiate = [=]() {
//        RayPlaneParams P;
//        P.p0 = { double(randi(-2, 2)), double(randi(-2, 2)), double(randi(0, 4)) };
//        P.n = { double(randi(0, 1)), double(randi(1, 2)), double(randi(1, 2)) };
//        P.d = { double(randi(-2,2)), double(randi(-2,2)), double(randi(-1,2)) };
//        if (std::fabs(P.d[0]) + std::fabs(P.d[1]) + std::fabs(P.d[2]) < 1e-9) P.d = { 0,0,1 };
//        P.o = { double(randi(-2, 2)), double(randi(-2, 2)), double(randi(-2, 2)) };
//
//        auto R = intersect_pointnormal(P.o, P.d, P.p0, P.n, EPS);
//
//        ostringstream q, a; q.setf(std::ios::fixed); a.setf(std::ios::fixed);
//        q << "Ray r(t)=o+t d, o=" << vec3s(P.o) << ", d=" << vec3s(P.d)
//            << "; Plane (x - p0)·n = 0 with p0=" << vec3s(P.p0) << ", n=" << vec3s(P.n) << ".\n"
//            << "Tasks: (1) Classify; (2) If hit, compute nearest non-negative t and intersection; "
//            << "(3) Return the **unit** plane normal (consistent with n). ε=" << EPS << ".";
//
//        a << setprecision(6);
//        if (R.state == HitState::Parallel) a << "Parallel.";
//        else if (R.state == HitState::Coplanar) a << "Coplanar.";
//        else if (R.state == HitState::Behind) a << "t=" << R.t << " < 0 (behind origin). unit normal=" << vec3s(R.unitN, 6);
//        else if (R.state == HitState::Hit) a << "t=" << R.t << ", hit=" << vec3s(R.hit, 3) << ", unit normal=" << vec3s(R.unitN, 6);
//        else a << "No hit.";
//        return Generated{ q.str(), a.str(), mkMetaJSON(T.name, T.rubric, Difficulty::D4, P) };
//        };
//    return T;
//}
//
//// ---------------- Bank & selector ----------------
//static vector<RayPlaneTemplate> buildRayPlaneBank() {
//    vector<RayPlaneTemplate> v;
//    v.push_back(T_D1_axisZ());
//    v.push_back(T_D2_general_ndotx_plus_d());
//    v.push_back(T_D3_pointNormal_classify_nearestT());
//    v.push_back(T_D4_pointNormal_with_unitNormal());
//    return v;
//}
//static Generated generateRayPlane(Difficulty D) {
//    auto bank = buildRayPlaneBank();
//    auto band = bandFor(D);
//    vector<RayPlaneTemplate> pool;
//    for (auto& t : bank) {
//        int s = t.rubric.total();
//        if (s >= band.first && s <= band.second) pool.push_back(t);
//    }
//    if (pool.empty()) pool = bank; // 极端回退
//    auto& T = pool[randi(0, int(pool.size()) - 1)];
//    return T.instantiate();
//}
//
//// ---------------- Interactive UI + export ----------------
//static int readInt(const string& prompt, int lo, int hi) {
//    while (true) {
//        cout << prompt;
//        int x; if (cin >> x && x >= lo && x <= hi) return x;
//        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
//        cout << "Invalid. Enter [" << lo << "," << hi << "].\n";
//    }
//}
//static unsigned readSeedOrRandom() {
//    cout << "Seed (enter -1 for random): ";
//    long long s; if (!(cin >> s)) { cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); s = -1; }
//    if (s < 0) { s = (long long)std::time(nullptr); }
//    return static_cast<unsigned>(s);
//}
//
//int main() {
//    cout << "=== Ray–Plane Generator (Understanding-first edition) ===\n\n";
//    bool running = true;
//    while (running) {
//        cout << "Menu\n 1) Generate Ray–Plane questions\n 0) Exit\n";
//        int op = readInt("Select: ", 0, 1);
//        if (op == 0) break;
//
//        unsigned seed = readSeedOrRandom(); reseed(seed);
//        int d = readInt("Difficulty 1-4: ", 1, 4);
//        Difficulty D = static_cast<Difficulty>(d);
//        int cnt = readInt("How many questions (1-20): ", 1, 20);
//
//        cout << "Export to output.jsonl ? 1=Yes 0=No : ";
//        int ex = 0; if (!(cin >> ex)) { cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); ex = 0; }
//        ofstream ofs; if (ex == 1) ofs.open("output.jsonl", ios::app);
//
//        cout << "\n--- Generated (seed=" << seed << ") ---\n";
//        for (int i = 1; i <= cnt; ++i) {
//            auto g = generateRayPlane(D);
//            cout << i << ". " << g.prompt << "\n\n";
//            cout << "   Ref: " << g.answer << "\n";
//            cout << "   Meta: " << g.metadata << "\n\n";
//            if (ofs) ofs << g.metadata << "\n"; // 一行一题
//        }
//        cout << "-------------------------------\n\n";
//    }
//    cout << "Bye.\n";
//    return 0;
//}
