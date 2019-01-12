#pragma once

// MIT License
//
// Copyright 2018 Michel Donais
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this
// software and associated documentation files (the "Software"), to deal in the Software
// without restriction, including without limitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be included in all copies
// or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// Attribution is given where appropriate. The only verbatim code parts are also under
// the MIT License:
// - https://github.com/fontello/cubic2quad

//
// The goal of this algorithm is to quickly convert a Cubic bezier to Quads.
//
// Sharing back the love.
//

#include <array>
#include <vector>

namespace CubicImpl
{
    // See https://stackoverflow.com/a/17728525/1060383
    inline constexpr int64_t ipow(int64_t base, int exp, int64_t result = 1) {
        return exp < 1 ? result : ipow(base*base, exp/2, (exp % 2) ? result*base : result);
    }
    
    // Simplest, no fuss, no verification vector implementation. Assumes you reserve first and never go out of bounds
    template <typename T, std::size_t capacity>
    struct stack_vector_impl
    {
        class allocator : public std::allocator<T>
        {
        public:
            using pointer = typename std::allocator<T>::pointer;
            using size_type = typename std::allocator<T>::size_type;
            
            using fake_array = typename std::array<T, capacity>;
            std::array<char, sizeof(fake_array)> stack_data;
            
            pointer allocate(size_type toAllocate, void* = 0)
            {
                return reinterpret_cast<pointer>(&stack_data);
            }
            void deallocate(pointer, size_type) {
            }
            
            // Copy constructor does nothing. Data will be copied by vector
            allocator() {}
            allocator(const allocator& rhs) {}
        };
    };
    template <typename T, std::size_t capacity>
    using stack_vector = std::vector<T, typename stack_vector_impl<T, capacity>::allocator>;
}

template <typename T, typename InnerFloat = float>
struct Cubic
{
    struct Point
    {
        T x, y;
        
        constexpr Point operator + (const Point& p) const
        {
            return { x + p.x, y + p.y };
        }
        constexpr Point operator - (const Point& p) const
        {
            return { x - p.x, y - p.y };
        }
        template <typename Oper>
        constexpr Point operator * (Oper f) const
        {
            return { x * f, y * f };
        }
        constexpr InnerFloat operator * (const Point& p) const
        {
            return x*p.x + y*p.y;
        }
        template <typename Oper>
        constexpr Point operator / (Oper f) const
        {
            return { x / f, y / f };
        }
        constexpr InnerFloat square() const
        {
            return x*x + y*y;
        }
        constexpr InnerFloat dist() const
        {
            return sqrtf(square());
        }
        constexpr Point min(const Point& p) const
        {
            return { std::min(x, p.x), std::min(y, p.y) };
        }
        constexpr Point max(const Point& p) const
        {
            return { std::max(x, p.x), std::max(y, p.y) };
        }
    };
    struct PointQuad {
        Point p,q;
    };
    
    // Bigger this value is, less number of points we will have, but worst the precision will be
    constexpr static const InnerFloat rangeEpsilon = std::numeric_limits<InnerFloat>::is_integer ? (InnerFloat)1 : (InnerFloat)(std::numeric_limits<InnerFloat>::epsilon() * 100000);
    
    // Every new cubic-to-quad segment will be verified at 4 different locations for adequate distance: 2 at left, 2 at right.
    constexpr static const auto numVerifyPointsPerCubic = 4;
    
    // It's possible to get 3 inflections, and with this value, every inflection can get 4 points, giving a maximum 12 quads, or 25 points total.
    constexpr static const auto maxRecurseSplitLevel = 2;
    
    using PointQuadVec = CubicImpl::stack_vector<PointQuad, 6 * CubicImpl::ipow(2, maxRecurseSplitLevel) + 1>;
    
    Point p1, c1, c2, p2;
    
    
    // The algorithm used for CubicToQuads is geared towards making calculations as fast as possible.
    // In order to make sure the quads are adequately placed, system checks at a few different places whether we are conforming
    // to the curve.
    
    // Returns a vector of points with quads, where the last one is always p=q.
    PointQuadVec toQuad() const
    {
        // Initial split by inflection
        auto splitCubic(getSplitCubicByInflection());
        
        const auto size = (p1.max(c1.max(c2.max(p2))) - p1.min(c1.min(c2.min(p2)))).dist();
        const auto epsilon = rangeEpsilon * splitCubic.size() * size * 100;
        
        // Further subdivisions where required, otherwise copy
        PointQuadVec result;
        for (auto& cubic : splitCubic)
        {
            recurseRange(result, maxRecurseSplitLevel, epsilon, cubic);
        }
        result.push_back({splitCubic.back().p2, splitCubic.back().p2});
        return result;
    }
    
    
protected:
    struct CubicQuad;
    
    using CubicVec = CubicImpl::stack_vector<CubicQuad, 6>;
    using FloatVec = CubicImpl::stack_vector<InnerFloat, 5>;
    
    static void pushIfValidRange(FloatVec& vec, const InnerFloat& value)
    {
        if (value > rangeEpsilon && value < (1 - rangeEpsilon))
        {
            vec.push_back(value);
        }
    }
    constexpr static bool isZero(InnerFloat value)
    {
        return abs(value) < rangeEpsilon;
    }
    
    static void recurseRange(PointQuadVec& vec, int level, InnerFloat epsilon, CubicQuad& rhs)
    {
        if (level > 0 && rhs.dist > epsilon)
        {
            CubicQuad lhs = rhs.splitAt(0.5);
            recurseRange(vec, level - 1, epsilon * 2, lhs);
            rhs.update();
            recurseRange(vec, level - 1, epsilon * 2, rhs);
        }
        else
        {
            if (level == 0 && rhs.dist > epsilon)
            {
                //printf("Not enough range @%f of %f\n", rhs.dist, epsilon);
            }
            vec.push_back({rhs.p1, rhs.quadC});
        }
    }
    
    static void solveQuad(InnerFloat a, InnerFloat b, InnerFloat c, FloatVec& out)
    {
        // Return to high school and solve (-b+-sqrt(b2-4ac))/2a
        if (a == 0)
        {
            if (b != 0)
            {
                pushIfValidRange(out, -c/b);
            }
        }
        else
        {
            const auto a2 = a * 2;
            const auto d = b*b - 4*a*c;
            
            if (isZero(d))
            {
                pushIfValidRange(out, -b / a2);
            }
            else if (d > 0)
            {
                const auto dSqrt = sqrtf(d);
                if (a < 0.f)            // Sort from smallest to biggest
                {
                    pushIfValidRange(out, (-b + dSqrt) / a2);
                    pushIfValidRange(out, (-b - dSqrt) / a2);
                }
                else
                {
                    pushIfValidRange(out, (-b - dSqrt) / a2);
                    pushIfValidRange(out, (-b + dSqrt) / a2);
                }
            }
        }
    }
    
    // This function and getMaxDist are pretty much copied from https://github.com/fontello/cubic2quad
    static void solveCubic(InnerFloat a, InnerFloat b, InnerFloat c, InnerFloat d, FloatVec& out)
    {
        // Return to college and solve our closest pal to a*x^3 + b*x^2 + c*x + d = 0
        if (isZero(a))
        {
            return solveQuad(b, c, d, out);
        }
        
        const InnerFloat xn = -b / (3*a); // point of symmetry x coordinate
        const InnerFloat yn = ((a * xn + b) * xn + c) * xn + d; // point of symmetry y coordinate
        const InnerFloat deltaSq = (b*b - 3*a*c) / (9*a*a); // delta^2
        const InnerFloat hSq = 4*a*a * deltaSq*deltaSq*deltaSq; // h^2
        const InnerFloat D3 = yn*yn - hSq;
        
        if (isZero(D3))
        { // 2 real roots
            const auto delta1 = cbrt(yn/(2*a));
            pushIfValidRange(out, xn - 2 * delta1);
            pushIfValidRange(out, xn + delta1);
        }
        else if (D3 > 0)
        { // 1 real root
            const auto D3Sqrt = sqrt(D3);
            pushIfValidRange(out, xn + cbrt((-yn + D3Sqrt)/(2*a)) + cbrt((-yn - D3Sqrt)/(2*a)));
        }
        else
        { // 3 real roots
            const InnerFloat theta = acos(-yn / sqrt(hSq)) / 3;
            const InnerFloat delta = sqrt(deltaSq);
            pushIfValidRange(out, xn + 2 * delta * cos(theta));
            pushIfValidRange(out, xn + 2 * delta * cos(theta + (InnerFloat)M_PI * 2 / 3));
            pushIfValidRange(out, xn + 2 * delta * cos(theta + (InnerFloat)M_PI * 4 / 3));
        }
    }
    
    struct PowerCoefficients
    {
        Point a, b, c, d;
        
        constexpr Point solve(InnerFloat t) const
        {
            // Solve a*t^3 + b*t^2 + c*t + d
            return ((a*t + b)*t + c)*t + d;
        }
        constexpr Point derivative(InnerFloat t) const
        {
            // Solve d/dt(solve(t))
            return (a*t*3 + b*2)*t + c;
        }
    };
    constexpr PowerCoefficients toPowerCoefficients() const
    {
        return
        {
            p2 - p1 + (c1 - c2) * 3,
            (p1 + c2) * 3 - c1 * 6,
            (c1 - p1) * 3,
            p1
        };
    }
    
    // Based on http://www.caffeineowl.com/graphics/2d/vectorial/cubic-inflexion.html.
    // Goal is to solve (-b+-sqrt(b2-4ac))/2a
    constexpr InnerFloat getA() const
    {
        return
          p1.x * (                c1.y - 2 * c2.y +     p2.y)
        - c1.x * (     p1.y            - 3 * c2.y + 2 * p2.y)
        + c2.x * ( 2 * p1.y - 3 * c1.y            +     p2.y)
        - p2.y * (     p1.y - 2 * c1.y +     c2.y           );
    }
    constexpr InnerFloat getB() const
    {
        return
        - p1.x * (            2 * c1.y - 3 * c2.y +     p2.y)
        + c1.x * ( 2 * p1.y            - 3 * c2.y +     p2.y)
        - c2.x * ( 3 * p1.y - 3 * c1.y                      )
        + p2.x * (     p1.y -     c1.y                      );
    }
    constexpr InnerFloat getC() const
    {
        return
          p1.x * (                c1.y -     c2.y           )
        + c1.x * (-    p1.y            +     c2.y           )
        + c2.x * (     p1.y -     c1.y                      );
    }
    
    // This is a glorified cubic with precalculated single-quad data
    struct CubicQuad : public Cubic<T, InnerFloat>
    {
        using Cubic = Cubic<T, InnerFloat>;
        
        PowerCoefficients coef;
        Point f1;
        Point f2;
        Point f1d;
        Point f2d;
        InnerFloat d;
        
        bool isStraightLine;
        Point quadC;
        float dist;
        
        CubicQuad(const Cubic& c):
        Cubic(c)
        {
            update();
        }
        
        CubicQuad(const Cubic&& c):
        Cubic(c)
        {
            update();
        }
        
        CubicQuad(const Point& p1_, const Point& c1_, const Point& c2_, const Point& p2_):
        Cubic({p1_, c1_, c2_, p2_})
        {
            update();
        }
        
        void update()
        {
            coef = toPowerCoefficients();
            f1 = coef.solve(0);
            f2 = coef.solve(1);
            f1d = coef.derivative(0);
            f2d = coef.derivative(1);
            d = segmentDenominator();
            isStraightLine = isZero(d);
            quadC = getQuadC();
            dist = getMaxDist();
        }
        
    protected:
        constexpr InnerFloat segmentDenominator() const
        {
            return f2d.x*f1d.y - f1d.x*f2d.y;
        }
        
        constexpr Point straightLineSegment() const
        {
            return (f1 + f2) / 2.f;
        }
        
    private:
        constexpr InnerFloat curvedLineSegmentE1() const
        {
            return f2.y*f2d.x - f2.x*f2d.y;
        }
        constexpr InnerFloat curvedLineSegmentE2() const
        {
            return f1.x*f1d.y - f1.y*f1d.x;
        }
    protected:
        constexpr Point curvedLineSegment() const
        {
            return Point({
                f1d.x*(curvedLineSegmentE1()) + f2d.x*(curvedLineSegmentE2()),
                f1d.y*(curvedLineSegmentE1()) + f2d.y*(curvedLineSegmentE2())
            }) / d;
        }
        
        constexpr Point getQuadC() const
        {
            return isStraightLine ? straightLineSegment() : curvedLineSegment();
        }
        
        // This function and solveCubic are pretty much copied from https://github.com/fontello/cubic2quad
        InnerFloat getMaxDist() const
        {
            InnerFloat maxDist = 0;
            
            constexpr const InnerFloat dt(1.f/(numVerifyPointsPerCubic + 1.f));
            InnerFloat t = dt;
            for (int i=0; i < numVerifyPointsPerCubic; ++i)
            {
                const auto cubicpt = coef.solve(t);
                t += dt;
                
                const auto a = p1 + p2 - quadC * 2;
                const auto b = (quadC - p1) * 2;
                const auto c = p1;
                const auto cSubPt = c - cubicpt;
                const auto e3 = a.square() * 2;
                const auto e2 = a * b * 3;
                const auto e1 = b.square() + a * cSubPt * 2;
                const auto e0 = cSubPt * b;
                
                FloatVec candidates;
                solveCubic(e3, e2, e1, e0, candidates);
                candidates.push_back(0);
                candidates.push_back(1);
                
                InnerFloat minDist = MAXFLOAT;
                for (auto candidateT : candidates)
                {
                    const auto distance = ((a*candidateT + b)*candidateT + c - cubicpt).dist();
                    if (distance < minDist)
                    {
                        minDist = distance;
                        if (distance == 0) break;
                    }
                }
                
                if (minDist > maxDist)
                {
                    maxDist = minDist;
                }
            }
            return maxDist;
        }
    };
    
    Cubic splitAt(float pos)
    {
        Cubic leftSide;
        
        const auto u = pos;
        const auto v = 1-u;
        
        leftSide.p1   = p1;
        leftSide.c1   = p1 * u          + c1 * v;
        const Point s = c1 * u          + c2 * v;
        leftSide.c2   = leftSide.c1 * u + s  * v;
        
        c2            = c2 * u          + p2 * v;
        c1            = s * u           + c2 * v;
        p1            = leftSide.c2 * u + c1 * v;
        leftSide.p2   = p1;
        
        return leftSide;
    }
    
    CubicVec getSplitCubicByInflection() const
    {
        const auto a = getA();
        const auto b = getB();
        const auto c = getC();
        FloatVec inflectionPosVec;
        solveQuad(a, b, c, inflectionPosVec);
        
        CubicVec result;
        result.reserve(inflectionPosVec.size() + 1);
        
        if (inflectionPosVec.empty())
        {
            result.push_back(std::move(CubicQuad(*this)));
            return result;
        }
        
        // Split at inflectionPos https://math.stackexchange.com/questions/877725
        Cubic rightSide(*this);
        float prevPoint = 0;
        for (auto inflectionPos : inflectionPosVec)
        {
            float pos = (1 - inflectionPos) / (1 - prevPoint);
            result.push_back(std::move(CubicQuad(std::move(rightSide.splitAt(pos)))));
            prevPoint = inflectionPos;
        }
        result.push_back(std::move(CubicQuad(std::move(rightSide))));
        
        return result;
    }
};

