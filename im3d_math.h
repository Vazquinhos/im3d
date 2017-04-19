#pragma once
#ifndef im3d_math_h
#define im3d_math_h

// im3d_math.h is optional - include only if you want to use the Im3d math types directly

#include "im3d.h"

#include <cmath>

namespace Im3d
{
    // Common math operators
    template <typename T> inline T Max(T _a, T _b) { return _a < _b ? _b : _a; }
    template <typename T> inline T Min(T _a, T _b) { return _a < _b ? _a : _b; }
    template <typename T> inline T Clamp(T _a, T _min, T _max) { return Max(Min(_a, _max), _min); }
    template <typename T> inline T Sqrt(T _val) { return std::sqrt(_val); }
    template <typename T> inline T Sin(T _ang) { return std::sin(_ang); }
    template <typename T> inline T Cos(T _ang) { return std::cos(_ang); }
    template <typename T> inline T Tan(T _ang) { return std::tan(_ang); }
    template <typename T> inline T SinCos(T _ang, T& s_, T& c_) { s_ = Sin(_ang);  c_ = Cos(_ang); }
    template <typename T> inline T ASin(T _val) { return std::asin(_val); }
    template <typename T> inline T ACos(T _val) { return std::acos(_val); }
    template <typename T> inline T ATan(T _val) { return std::atan(_val); }
    template <typename T> inline T ATan2(T _x, T _y) { return std::atan2(_x, _y); }
    template <typename T> inline T Square(T _val) { return _val*_val; }
    template <typename T> inline T PowN(T _val, int _n) { T r = 1.0f; if (_n > 0) for (int i = 0; i < _n; ++i) { r *= _val; } else if (_n < 0) { for (int i = 0; i < -_n; ++i) r /= _val; } return r; }
    template <typename T> inline T Log(T _val) { return std::log(_val); }
    template <typename T> inline T LogBase(T _val, T base) { return T(Log(_val) / Log(base)); }
    template <typename T> inline T Abs(T _val) { return std::fabs(_val); }
    template <typename T> inline T Floor(T _val) { return std::floor(_val); }
    template <typename T> inline T Ceil(T _val) { return std::ceil(_val); }
    template <typename T> inline T Round(T _val) { return (_val > 0.0f) ? std::floorf(_val + 0.5f) : std::ceilf(_val - 0.5f); }
    template <typename T> inline T Exp(T _val) { return std::exp(_val); }
    template <typename T> inline T Lerp(const T _x, T _y, T _s) { return _x + (_y - _x) * _s; }
    template <typename T> inline T Sign(T _val) { int i = static_cast<int>((*(int*)(&_val) & 0x80000000 | 0x3F800000)); return *(T*)(&i); }
    template <typename T> inline T Saturate(T _val) { return Clamp(_val, 0.0f, 1.0f); }
    template <typename T> inline T Smoothstep(T start, T end, T x) { x = Saturate((x - start) / (end - start)); return x * x * (3.0f - 2.0f * x); }
    template <typename T> inline bool IsPositive(T _val) { return (*(int*)(&_val) & 0x80000000) == 0; }

    // Vec2
    inline Vec2  operator+(const Vec2& _lhs, const Vec2& _rhs) { return Vec2(_lhs.x + _rhs.x, _lhs.y + _rhs.y); }
    inline Vec2  operator-(const Vec2& _lhs, const Vec2& _rhs) { return Vec2(_lhs.x - _rhs.x, _lhs.y - _rhs.y); }
    inline Vec2  operator*(const Vec2& _lhs, const Vec2& _rhs) { return Vec2(_lhs.x * _rhs.x, _lhs.y * _rhs.y); }
    inline Vec2  operator/(const Vec2& _lhs, const Vec2& _rhs) { return Vec2(_lhs.x / _rhs.x, _lhs.y / _rhs.y); }
    inline Vec2  operator*(const Vec2& _lhs, float _rhs) { return Vec2(_lhs.x * _rhs, _lhs.y * _rhs); }
    inline Vec2  operator/(const Vec2& _lhs, float _rhs) { return Vec2(_lhs.x / _rhs, _lhs.y / _rhs); }
    inline Vec2  operator-(const Vec2& _v) { return Vec2(-_v.x, -_v.y); }
    inline float Dot(const Vec2& _lhs, const Vec2& _rhs) { return _lhs.x * _rhs.x + _lhs.y * _rhs.y; }
    inline float Length(const Vec2& _v) { return sqrtf(Dot(_v, _v)); }
    inline float Length2(const Vec2& _v) { return Dot(_v, _v); }
    inline Vec2  Abs(const Vec2& _v) { return Vec2(fabs(_v.x), fabs(_v.y)); }
    inline Vec2  Normalize(const Vec2& _v) { return _v / Length(_v); }
    inline Vec2 Slerp(const Vec2& start, const Vec2& end, float percent)
    {
        Vec2 lStart = Normalize(start);
        Vec2 lEnd = Normalize(end);

        float dot = Dot(lStart, lEnd);
        dot = Clamp(dot, -1.0f, 1.0f);
        float theta = acos(dot) * percent;
        Vec2 relative_vector = lEnd - lStart * dot;
        relative_vector = Normalize(relative_vector);
        return ((lStart * Cos(theta)) + (relative_vector * Sin(theta)));
    }
    inline Vec2 Lerp(const Vec2& start, const Vec2& end, float percent)
    {
        return Vec2( Lerp(start.x, end.x, percent), Lerp(start.y, end.y, percent) );
    }
    inline Vec2 Perpendicular(const Vec2& v)
    {
        return Vec2(v.y, -v.x);
    }
    static inline Vec2 Reflect(const Vec2& d, const Vec2& n)
    {
        Vec2 norm_n = Normalize(n);
        return d - norm_n * Dot(norm_n, d) * 2.0f;
    }
    inline float Angle(const Vec2& a, const Vec2& b)
    {
        return ATan2(a.x* b.y - a.y*b.x, Dot(a, b));
    }

    // Vec3
    inline Vec3  operator+(const Vec3& _lhs, const Vec3& _rhs) { return Vec3(_lhs.x + _rhs.x, _lhs.y + _rhs.y, _lhs.z + _rhs.z); }
    inline Vec3  operator-(const Vec3& _lhs, const Vec3& _rhs) { return Vec3(_lhs.x - _rhs.x, _lhs.y - _rhs.y, _lhs.z - _rhs.z); }
    inline Vec3  operator*(const Vec3& _lhs, const Vec3& _rhs) { return Vec3(_lhs.x * _rhs.x, _lhs.y * _rhs.y, _lhs.z * _rhs.z); }
    inline Vec3  operator/(const Vec3& _lhs, const Vec3& _rhs) { return Vec3(_lhs.x / _rhs.x, _lhs.y / _rhs.y, _lhs.z / _rhs.z); }
    inline Vec3  operator*(const Vec3& _lhs, float _rhs) { return Vec3(_lhs.x * _rhs, _lhs.y * _rhs, _lhs.z * _rhs); }
    inline Vec3  operator/(const Vec3& _lhs, float _rhs) { return Vec3(_lhs.x / _rhs, _lhs.y / _rhs, _lhs.z / _rhs); }
    inline Vec3  operator-(const Vec3& _v) { return Vec3(-_v.x, -_v.y, -_v.z); }
    inline float Dot(const Vec3& _lhs, const Vec3& _rhs) { return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z; }
    inline float Length(const Vec3& _v) { return sqrtf(Dot(_v, _v)); }
    inline float Length2(const Vec3& _v) { return Dot(_v, _v); }
    inline Vec3  Abs(const Vec3& _v) { return Vec3(fabs(_v.x), fabs(_v.y), fabs(_v.z)); }
    inline Vec3  Normalize(const Vec3& _v) { return _v / Length(_v); }
    inline Vec3  Cross(const Vec3& _a, const Vec3& _b)
    {
        return Vec3(
            _a.y * _b.z - _b.y * _a.z,
            _a.z * _b.x - _b.z * _a.x,
            _a.x * _b.y - _b.x * _a.y
        );
    }
    inline Vec3 Slerp(const Vec3& start, const Vec3& end, float percent)
    {
        Vec3 lStart = Normalize(start);
        Vec3 lEnd = Normalize(end);

        float dot = Dot(lStart, lEnd);
        dot = Clamp(dot, -1.0f, 1.0f);
        float theta = acos(dot) * percent;
        Vec3 relative_vector = lEnd - lStart * dot;
        relative_vector = Normalize(relative_vector);
        return ((lStart * cos(theta)) + (relative_vector * sin(theta)));
    }
    inline Vec3 Lerp(const Vec3& start, const Vec3& end, float percent)
    {
        return Vec3(Lerp(start.x, end.x, percent), Lerp(start.y, end.y, percent), Lerp(start.z, end.z, percent));
    }
    inline Vec3 Perpendicular(const Vec3& v)
    {
        Vec3 perp;

        float x = Abs(v.x);
        float y = Abs(v.y);
        float z = Abs(v.z);
        float minVal = Min(x, y);
        minVal = Min(minVal, z);

        if (minVal == x)      perp = Cross(v, Vec3::right);
        else if (minVal == y) perp = Cross(v, Vec3::up);
        else                  perp = Cross(v, Vec3::forward);

        return Normalize(perp);
    }
    static inline Vec3 ProjectOnPlane(const Vec3& vector, const Vec3& planeNormal)
    {
        return vector + planeNormal * -Dot(planeNormal, vector);
    }
    static inline Vec3 Reflect(const Vec3& d, const Vec3& n)
    {
        Vec3 norm_n = Normalize(n);
        return d - norm_n * Dot(norm_n, d) * 2.0f;
    }
    inline float Angle(const Vec3& a, const Vec3& b)
    {
        return ATan2(Length2( Cross( a, b ) ), Dot(a, b));
    }

    // Vec4
    inline Vec4  operator+(const Vec4& _lhs, const Vec4& _rhs) { return Vec4(_lhs.x + _rhs.x, _lhs.y + _rhs.y, _lhs.z + _rhs.z, _lhs.w + _rhs.w); }
    inline Vec4  operator-(const Vec4& _lhs, const Vec4& _rhs) { return Vec4(_lhs.x - _rhs.x, _lhs.y - _rhs.y, _lhs.z - _rhs.z, _lhs.w - _rhs.w); }
    inline Vec4  operator*(const Vec4& _lhs, const Vec4& _rhs) { return Vec4(_lhs.x * _rhs.x, _lhs.y * _rhs.y, _lhs.z * _rhs.z, _lhs.w * _rhs.w); }
    inline Vec4  operator/(const Vec4& _lhs, const Vec4& _rhs) { return Vec4(_lhs.x / _rhs.x, _lhs.y / _rhs.y, _lhs.z / _rhs.z, _lhs.w / _rhs.w); }
    inline Vec4  operator*(const Vec4& _lhs, float _rhs) { return Vec4(_lhs.x * _rhs, _lhs.y * _rhs, _lhs.z * _rhs, _lhs.w * _rhs); }
    inline Vec4  operator/(const Vec4& _lhs, float _rhs) { return Vec4(_lhs.x / _rhs, _lhs.y / _rhs, _lhs.z / _rhs, _lhs.w / _rhs); }
    inline Vec4  operator-(const Vec4& _v) { return Vec4(-_v.x, -_v.y, -_v.z, -_v.w); }
    inline float Dot(const Vec4& _lhs, const Vec4& _rhs) { return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z + _lhs.w * _rhs.w; }
    inline float Length(const Vec4& _v) { return sqrtf(Dot(_v, _v)); }
    inline float Length2(const Vec4& _v) { return Dot(_v, _v); }
    inline Vec4  Abs(const Vec4& _v) { return Vec4(fabs(_v.x), fabs(_v.y), fabs(_v.z), fabs(_v.w)); }
    inline Vec4  Normalize(const Vec4& _v) { return _v / Length(_v); }

    // Mat3
    inline Mat3 operator*(const Mat3& _lhs, const Mat3& _rhs)
    {
        Mat3 ret;
        ret(0, 0) = _lhs(0, 0) * _rhs(0, 0) + _lhs(0, 1) * _rhs(1, 0) + _lhs(0, 2) * _rhs(2, 0);
        ret(0, 1) = _lhs(0, 0) * _rhs(0, 1) + _lhs(0, 1) * _rhs(1, 1) + _lhs(0, 2) * _rhs(2, 1);
        ret(0, 2) = _lhs(0, 0) * _rhs(0, 2) + _lhs(0, 1) * _rhs(1, 2) + _lhs(0, 2) * _rhs(2, 2);
        ret(1, 0) = _lhs(1, 0) * _rhs(0, 0) + _lhs(1, 1) * _rhs(1, 0) + _lhs(1, 2) * _rhs(2, 0);
        ret(1, 1) = _lhs(1, 0) * _rhs(0, 1) + _lhs(1, 1) * _rhs(1, 1) + _lhs(1, 2) * _rhs(2, 1);
        ret(1, 2) = _lhs(1, 0) * _rhs(0, 2) + _lhs(1, 1) * _rhs(1, 2) + _lhs(1, 2) * _rhs(2, 2);
        ret(2, 0) = _lhs(2, 0) * _rhs(0, 0) + _lhs(2, 1) * _rhs(1, 0) + _lhs(2, 2) * _rhs(2, 0);
        ret(2, 1) = _lhs(2, 0) * _rhs(0, 1) + _lhs(2, 1) * _rhs(1, 1) + _lhs(2, 2) * _rhs(2, 1);
        ret(2, 2) = _lhs(2, 0) * _rhs(0, 2) + _lhs(2, 1) * _rhs(1, 2) + _lhs(2, 2) * _rhs(2, 2);
        return ret;
    }
    inline Vec3 operator*(const Mat3& _m, const Vec3& _v)
    {
        return Vec3(
            _m(0, 0) * _v.x + _m(0, 1) * _v.y + _m(0, 2) * _v.z,
            _m(1, 0) * _v.x + _m(1, 1) * _v.y + _m(1, 2) * _v.z,
            _m(2, 0) * _v.x + _m(2, 1) * _v.y + _m(2, 2) * _v.z
        );
    }
    inline Vec4 operator*(const Mat3& _m, const Vec4& _v)
    {
        return Vec4(
            _m(0, 0) * _v.x + _m(0, 1) * _v.y + _m(0, 2) * _v.z,
            _m(1, 0) * _v.x + _m(1, 1) * _v.y + _m(1, 2) * _v.z,
            _m(2, 0) * _v.x + _m(2, 1) * _v.y + _m(2, 2) * _v.z,
            _v.w
        );
    }
    Mat3 Transpose(const Mat3& _m);
    Vec3 ToEulerXYZ(const Mat3& _m);
    Mat3 FromEulerXYZ(Vec3& _xyz);
    Mat3 Rotation(const Vec3& _axis, float _rads); // _axis must be unit length
    Mat3 Scale(const Vec3& _s);

    // Mat4
    inline Mat4 operator*(const Mat4& _lhs, const Mat4& _rhs)
    {
        Mat4 ret;
        ret(0, 0) = _lhs(0, 0) * _rhs(0, 0) + _lhs(0, 1) * _rhs(1, 0) + _lhs(0, 2) * _rhs(2, 0) + _lhs(0, 3) * _rhs(3, 0);
        ret(0, 1) = _lhs(0, 0) * _rhs(0, 1) + _lhs(0, 1) * _rhs(1, 1) + _lhs(0, 2) * _rhs(2, 1) + _lhs(0, 3) * _rhs(3, 1);
        ret(0, 2) = _lhs(0, 0) * _rhs(0, 2) + _lhs(0, 1) * _rhs(1, 2) + _lhs(0, 2) * _rhs(2, 2) + _lhs(0, 3) * _rhs(3, 2);
        ret(0, 3) = _lhs(0, 0) * _rhs(0, 3) + _lhs(0, 1) * _rhs(1, 3) + _lhs(0, 2) * _rhs(2, 3) + _lhs(0, 3) * _rhs(3, 3);
        ret(1, 0) = _lhs(1, 0) * _rhs(0, 0) + _lhs(1, 1) * _rhs(1, 0) + _lhs(1, 2) * _rhs(2, 0) + _lhs(1, 3) * _rhs(3, 0);
        ret(1, 1) = _lhs(1, 0) * _rhs(0, 1) + _lhs(1, 1) * _rhs(1, 1) + _lhs(1, 2) * _rhs(2, 1) + _lhs(1, 3) * _rhs(3, 1);
        ret(1, 2) = _lhs(1, 0) * _rhs(0, 2) + _lhs(1, 1) * _rhs(1, 2) + _lhs(1, 2) * _rhs(2, 2) + _lhs(1, 3) * _rhs(3, 2);
        ret(1, 3) = _lhs(1, 0) * _rhs(0, 3) + _lhs(1, 1) * _rhs(1, 3) + _lhs(1, 2) * _rhs(2, 3) + _lhs(1, 3) * _rhs(3, 3);
        ret(2, 0) = _lhs(2, 0) * _rhs(0, 0) + _lhs(2, 1) * _rhs(1, 0) + _lhs(2, 2) * _rhs(2, 0) + _lhs(2, 3) * _rhs(3, 0);
        ret(2, 1) = _lhs(2, 0) * _rhs(0, 1) + _lhs(2, 1) * _rhs(1, 1) + _lhs(2, 2) * _rhs(2, 1) + _lhs(2, 3) * _rhs(3, 1);
        ret(2, 2) = _lhs(2, 0) * _rhs(0, 2) + _lhs(2, 1) * _rhs(1, 2) + _lhs(2, 2) * _rhs(2, 2) + _lhs(2, 3) * _rhs(3, 2);
        ret(2, 3) = _lhs(2, 0) * _rhs(0, 3) + _lhs(2, 1) * _rhs(1, 3) + _lhs(2, 2) * _rhs(2, 3) + _lhs(2, 3) * _rhs(3, 3);
        ret(3, 0) = _lhs(3, 0) * _rhs(0, 0) + _lhs(3, 1) * _rhs(1, 0) + _lhs(3, 2) * _rhs(2, 0) + _lhs(3, 3) * _rhs(3, 0);
        ret(3, 1) = _lhs(3, 0) * _rhs(0, 1) + _lhs(3, 1) * _rhs(1, 1) + _lhs(3, 2) * _rhs(2, 1) + _lhs(3, 3) * _rhs(3, 1);
        ret(3, 2) = _lhs(3, 0) * _rhs(0, 2) + _lhs(3, 1) * _rhs(1, 2) + _lhs(3, 2) * _rhs(2, 2) + _lhs(3, 3) * _rhs(3, 2);
        ret(3, 3) = _lhs(3, 0) * _rhs(0, 3) + _lhs(3, 1) * _rhs(1, 3) + _lhs(3, 2) * _rhs(2, 3) + _lhs(3, 3) * _rhs(3, 3);
        return ret;
    }
    inline Vec3 operator*(const Mat4& _m, const Vec3& _pos)
    {
        return Vec3(
            _m(0, 0) * _pos.x + _m(0, 1) * _pos.y + _m(0, 2) * _pos.z + _m(0, 3),
            _m(1, 0) * _pos.x + _m(1, 1) * _pos.y + _m(1, 2) * _pos.z + _m(1, 3),
            _m(2, 0) * _pos.x + _m(2, 1) * _pos.y + _m(2, 2) * _pos.z + _m(2, 3)
        );
    }
    inline Vec4 operator*(const Mat4& _m, const Vec4& _v)
    {
        return Vec4(
            _m(0, 0) * _v.x + _m(0, 1) * _v.y + _m(0, 2) * _v.z + _m(0, 3) * _v.w,
            _m(1, 0) * _v.x + _m(1, 1) * _v.y + _m(1, 2) * _v.z + _m(1, 3) * _v.w,
            _m(2, 0) * _v.x + _m(2, 1) * _v.y + _m(2, 2) * _v.z + _m(2, 3) * _v.w,
            _m(3, 0) * _v.x + _m(3, 1) * _v.y + _m(3, 2) * _v.z + _m(3, 3) * _v.w
        );
    }
    Mat4 Inverse(const Mat4& _m);
    Mat4 Transpose(const Mat4& _m);
    Mat4 Translation(const Vec3& _t);
    Mat4 AlignZ(const Vec3& _axis, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f)); // generate an orthonormal bases with +z as _axis, which must be unit length
    Mat4 LookAt(const Vec3& _from, const Vec3& _to, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f)); // align _z with (_to - _from), set _from as translation

    struct Line
    {
        Vec3 m_origin;
        Vec3 m_direction; // unit length

        Line() {}
        Line(const Vec3& _origin, const Vec3& _direction);
    };
    struct Ray
    {
        Vec3 m_origin;
        Vec3 m_direction; // unit length

        Ray() {}
        Ray(const Vec3& _origin, const Vec3& _direction);
    };
    struct LineSegment
    {
        Vec3 m_start;
        Vec3 m_end;

        LineSegment() {}
        LineSegment(const Vec3& _start, const Vec3& _end);
    };
    struct Sphere
    {
        Vec3  m_origin;
        float m_radius;

        Sphere() {}
        Sphere(const Vec3& _origin, float _radius);
    };
    struct Plane
    {
        Vec3  m_normal;
        float m_offset;

        Plane() {}
        Plane(const Vec3& _normal, float _offset);
        Plane(const Vec3& _normal, const Vec3& _origin);
    };
    struct Capsule
    {
        Vec3  m_start;
        Vec3  m_end;
        float m_radius;

        Capsule() {}
        Capsule(const Vec3& _start, const Vec3& _end, float _radius);
    };


    // Ray-primitive intersections. Use Intersects() when you don't need t.
    bool  Intersects(const Ray& _ray, const Plane& _plane);
    bool  Intersect(const Ray& _ray, const Plane& _plane, float& t0_);
    bool  Intersects(const Ray& _ray, const Sphere& _sphere);
    bool  Intersect(const Ray& _ray, const Sphere& _sphere, float& t0_, float& t1_);
    bool  Intersects(const Ray& _ray, const Capsule& _capsule);
    bool  Intersect(const Ray& _ray, const Capsule& _capsule, float& t0_, float& t1_);

    // Find point t0_ along _line0 nearest to _line1 and t1_ along _line1 nearest to _line0.
    void  Nearest(const Line& _line0, const Line& _line1, float& t0_, float& t1_);
    // Find point tr_ along _ray nearest to _line and tl_ along _line nearest to _ray.
    void  Nearest(const Ray& _ray, const Line& _line, float& tr_, float& tl_);
    // Find point tr_ along _ray nearest to _segment, return point on segment nearest to _ray.
    Vec3  Nearest(const Ray& _ray, const LineSegment& _segment, float& tr_);

    float Distance2(const Ray& _ray, const LineSegment& _segment);

    extern const float Pi;
    extern const float TwoPi;
    extern const float HalfPi;

    inline float Radians(float _degrees) { return _degrees * Pi / 180.0f; }
    inline float Degrees(float _radians) { return _radians * 180.0f / Pi; }

    // Remap _x in [_start,_end] to [0,1].
    inline float Remap(float _x, float _start, float _end) { return Clamp(_x * (1.0f / (_end - _start)) + (-_start / (_end - _start)), 0.0f, 1.0f); }

    template <typename T>
    int Count();
    template <> inline int Count<Vec2>() { return 2; }
    template <> inline int Count<Vec3>() { return 3; }
    template <> inline int Count<Vec4>() { return 4; }
    template <> inline int Count<Mat4>() { return 16; }

    template <typename T>
    inline bool AllLess(const T& _a, const T& _b)
    {
        for (int i = 0, n = Count<T>(); i < n; ++i) {
            if (_a[i] > _b[i]) {
                return false;
            }
        }
        return true;
    }

} // namespace Im3d

#endif // im3d_math_h
