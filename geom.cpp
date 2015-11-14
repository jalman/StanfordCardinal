#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#define EPS 1e-9
#define PI 3.141592653589793
using namespace std;

// points are treated as vectors
struct point {
	double x, y;
	point(const point& other): x(other.x), y(other.y) {}
	point(double _x = 0, double _y = 0): x(_x), y(_y) {}
	point operator+(const point& other) const { return point(x + other.x, y + other.y); }
	point operator-(const point& other) const { return point(x - other.x, y - other.y); }
	point operator*(const double& other) const { return point(x * other, y * other); }
	point operator/(const double& other) const { return point(x / other, y / other); }
	point operator-() const { return point(-x, -y); }
	bool operator==(const point& other) const { return fabs(x - other.x) < EPS && fabs(y - other.y) < EPS; }
	bool operator!=(const point& other) const { return fabs(x - other.x) >= EPS || fabs(y - other.y) >= EPS; }
	inline double dot(const point& other) const { return x * other.x + y * other.y; }
	// z-component of cross product of (x1, y1, 0) and (x2, y2, 0), which gives signed area of paralellogram
	inline double cross(const point& other) const { return x * other.y - y * other.x; }
	inline double sqnorm() const { return x * x + y * y; } // squared norm = squared distance
	inline double norm() const { return sqrt(x * x + y * y); } // norm = length
	inline double angle() const { return atan2(y, x); } // angle in radians

	point normalize() const { // use this to normalize the point/vector
		double length = norm();
		return point(x / length, y / length);
	}
	point rotate(double angle) const { // rotate the point/vector
		double cs = cos(angle), sn = sin(angle);
		return point(cs * x - sn * y, sn * x + cs * y);
	}
	inline point perp() const { return point(-y, x); } // rotate by 90 degrees counterclockwise
};

// ccw < 0 - clockwise, ccw == 0 - collinear, ccw > 0 - counterclockwise. abs(ccw) = area of paralellogram.
inline double ccw(const point& a, const point& b, const point& c) { return (b - a).cross(c - b); }
inline bool collinear(const point& a, const point& b, const point& c) { return fabs(ccw(a, b, c)) < EPS; }

// lines are represented by line segments
struct line {
	point a, b;
	line(double _ax = 0, double _ay = 0, double _bx = 0, double _by = 0): a(_ax, _ay), b(_bx, _by) {}
	line(point _a, point _b): a(_a), b(_b) {}
	inline double length() const { return (b - a).norm(); }
	inline line reverse() const { return line(b, a); }
	inline double angle() const { return (b - a).angle(); } // angle of rotation of line
	// minor relative angle between two lines
	inline double angle(const line& other) const { return acos(fabs((b - a).normalize().dot((other.b - other.a).normalize()))); }
	inline point midpoint() const { return (a + b) * 0.5; }
	inline double projectionReference(const point& other) const { return (other - a).dot(b - a) / (b - a).sqnorm(); } // normalized coordinates
	point projection(const point& other) const { // projection of a point onto the line. projection = a + (b - a) * reference
		point dir = b - a;
		return a + dir * (other - a).dot(dir) / dir.sqnorm();
	}
	double subtendedAngle(const point& other) const { // angle subtended by line on point
		return acos((a - other).normalize().dot((b - other).normalize()));
	}
	point perpendicularFoot(const point& other) const { // perpendicular foot on line
		point dir = b - a;
		double lambda = dir.dot(other - a) / dir.sqnorm();
		return point(a.x + lambda * dir.x, a.y + lambda * dir.y);
	}
	double perpendicularDistance(const point& other) const { // perpendicular distance to line
		point foot = perpendicularFoot(other);
		return hypot(other.x - foot.x, other.y - foot.y);
	}
	point closestPoint(const point& other) const { // closest point in line segment
		point dir = b - a;
		double lambda = dir.dot(other - a) / dir.sqnorm();
		if (lambda >= 1 + EPS) return b;
		else if (lambda <= -EPS) return a;
		else return point(a.x + lambda * dir.x, a.y + lambda * dir.y);
	}
	double closestDistance(const point& other) const { // closest distance to line segment
		point closest = closestPoint(other);
		return hypot(other.x - closest.x, other.y - closest.y);
	}
	short pointPartition(const point& other) const { // 1 = left, 0 = middle, -1 = right
		double res = ccw(a, b, other);
		return (fabs(res) < EPS) ? 0 : (res > EPS ? 1 : -1);
	}
	inline bool parallel(const line& other) const {
		return collinear(point(), b - a, other.b - other.a);
	}
	inline bool collinearTo(const line& other) const { // are the two lines equivalent?
		return collinear(a, b, other.a) && collinear(a, b, other.b);
	}
	inline bool collinearTo(const point& other) const { // is the point on this line?
		return collinear(a, b, other);
	}
	bool intersects(const line& other) const { // true if the line segments intersect
		if (collinearTo(other)) {
			point dir = b - a;
			double length = dir.sqnorm(), mu = dir.dot(other.a - a) / length, lambda = dir.dot(other.b - a) / length;
			return (EPS < mu && mu + EPS < 1) || (EPS < lambda && lambda + EPS < 1) || (min(mu, lambda) <= EPS && max(mu, lambda) + EPS >= 1);
		} else return (ccw(a, b, other.a) * ccw(a, b, other.b) < -EPS) && (ccw(other.a, other.b, a) * ccw(other.a, other.b, b) < -EPS);
	}
	bool contacts(const line& other) const { // true if the line segments touch
		if (collinearTo(other)) {
			point dir = b - a;
			double length = dir.sqnorm(), mu = dir.dot(other.a - a) / length, lambda = dir.dot(other.b - a) / length;
			return (-EPS <= mu && mu < 1 + EPS) || (-EPS <= lambda && lambda < 1 + EPS) || (min(mu, lambda) <= EPS && max(mu, lambda) + EPS >= 1);
		} else return (ccw(a, b, other.a) * ccw(a, b, other.b) < EPS) && (ccw(other.a, other.b, a) * ccw(other.a, other.b, b) < EPS);
	}
	point intersectionPoint(const line& other) const { // gives intersection point of the two lines (extended beyond segments)
		if (parallel(other)) return point(); // the lines must not be parallel
		double abA = ccw(a, b, other.a), abB = ccw(a, b, other.b);
		if (fabs(abA - abB) < EPS) return a;
		double lambda = abA / (abA - abB);
		point dir = other.b - other.a;
		return point(other.a.x + dir.x * lambda, other.a.y + dir.y * lambda);
	}
};

// the circle is represented by its centre and its radius
struct circle {
	point c;
	double r;
	circle(point _c = point(), double _r = 0): c(_c), r(_r) {}
	circle(point _a, point _b, point _c) { // circle from three points on the circle
		c = point(), r = 0.0;
		if (!collinear(_a, _b, _c)) { // the three points must not be collinear
			point dab = (_b - _a).perp(), dbc = (_c - _b).perp();
			point mab = (_a + _b) * 0.5, mbc = (_b + _c) * 0.5;
			c = line(mab, mab + dab).intersectionPoint(line(mbc, mbc + dbc));
			r = (_a - c).norm();
		}
	}
	bool operator==(const circle& other) const { return c == other.c && fabs(r - other.r) < EPS; }
	bool operator!=(const circle& other) const { return c != other.c || fabs(r - other.r) >= EPS; }
	inline double area() const { return PI * r * r; }
	inline double circumference() const { return 2.0 * PI * r; }
	inline double sector(double angle) { return 0.5 * angle * r * r; }
	inline double segment(double angle) { return 0.5 * (angle - sin(angle)) * r * r; }
	inline double chord(double angle) { return 2.0 * r * r * (1.0 - cos(angle)); } // length of chord
	inline bool contains(const point& other) const { return (other - c).sqnorm() + EPS < r * r; }
	inline bool contacts(const point& other) const { return (other - c).sqnorm() < r * r + EPS; }
	inline bool contains(const line& other) const { return contains(other.a) && contains(other.b); }
	bool intersects(const line& other) const { return other.perpendicularDistance(c) < r + EPS; }
	bool contacts(const line& other) const { return other.perpendicularDistance(c) + EPS <= r; }
	int intersectsSegment(const line& other) const { // 0 - no, 1 - yes, 2 - segment cuts circle twice
		point foot = other.perpendicularFoot(c);
		double ref = other.projectionReference(foot);
		double dist = sqrt(r * r - (foot - c).sqnorm()) / other.length();
		return (ref + dist > EPS && ref + dist + EPS < 1) + (ref - dist > EPS && ref - dist + EPS < 1);
	}
	int contactsSegment(const line& other) const { // 0 - no, 1 - yes, 2 - segment touches circle twice
		point foot = other.perpendicularFoot(c);
		double ref = other.projectionReference(foot);
		double dist = sqrt(r * r - (foot - c).sqnorm()) / other.length();
		return (ref + dist > -EPS && ref + dist < 1 + EPS) + (ref - dist > -EPS && ref - dist < 1 + EPS);
	}
	point intersectionPoint(const line& other) const { // intersection between extended line and circle. reverse line to get other point.
		point foot = other.perpendicularFoot(c);
		point dir = foot - c;
		double dist = sqrt(r * r - dir.sqnorm());
		return foot + dir.perp().normalize() * dist;
	}
	inline bool contains(const circle& other) const { return (other.c - c).norm() + EPS < fabs(other.r - r); }
	inline bool intersects(const circle& other) const { return !contains(other) && (other.c - c).norm() + EPS < r + other.r; }
	inline bool contacts(const circle& other) const { return !contains(other) && (other.c - c).norm() < EPS + r + other.r; }
	point intersectionPoint(const circle& other) const { // the two circles must at least contact each other.
		if (!contacts(other)) return point();
		point dir = other.c - c;
		double a = (r * r - other.r * other.r + dir.sqnorm()) / (2.0 * dir.norm());
		double h = sqrt(r * r - a * a);
		point centre = c + dir.normalize() * a;
		return centre + dir.perp().normalize() * h;
	}
};

inline bool isChord(const line& l, double r) { // can this line be a chord of radius r circle?
	point dir = l.b - l.a;
	double determinant = r * r / dir.sqnorm() - 0.25;
	return determinant > -EPS;
}
circle circleFromChord(const line& l, double r) { // to get the other circle, reverse the line.
	point dir = l.b - l.a;
	double determinant = r * r / dir.sqnorm() - 0.25;
	double h = sqrt(determinant);
	return circle((l.a + l.b) * 0.5 + dir.perp() * h, r);
}

// polygons are lists of points with no repeated end vertices. use counterclockwise winding.
struct polygon {
	vector<point> p;
	polygon(vector<point> _p = vector<point>()): p(_p) {}
	inline point& operator[](const size_t& other) { return p[other]; }
	inline const point& operator[](const size_t& other) const { return p[other]; }
	inline bool empty() const { return p.empty(); }
	inline size_t size() const { return p.size(); }
	inline void push_back(const point& other) { p.push_back(other); }
	double perimeter() const { // perimeter of polygon
		double total = 0.0;
		for (size_t i = 0; i < p.size(); ++i)
			total += (p[(i + 1) % p.size()] - p[i]).norm();
		return total;
	}
	double area() const { // area of polygon
		double total = 0.0;
		for (size_t i = 0; i < p.size(); ++i)
			total += p[i].perp().dot(p[(i + 1) % p.size()]);
		return fabs(total) / 2.0;
	}
	bool isCCW() const { // is polygon a counterclockwise convex polygon?
		if (p.size() < 3) return false;
		for (size_t i = 0; i < p.size(); ++i)
			if (ccw(p[i], p[(i + 1) % p.size()], p[(i + 2) % p.size()]) < EPS) return false;
		return true;
	}
	bool isConvex() const { // is polygon convex
		if (p.size() < 3) return false;
		bool direction = ccw(p[0], p[1], p[2]) >= EPS;
		for (size_t i = 1; i < p.size(); ++i)
			if ((ccw(p[i], p[(i + 1) % p.size()], p[(i + 2) % p.size()]) >= EPS) != direction)
				return false;
		return true;
	}
	bool contains(const point& other) const { // contains point
		if (p.empty()) return false;
		double total = 0.0;
		for (size_t i = 0; i < p.size(); ++i)
			total += ((ccw(other, p[i], p[(i + 1) % p.size()]) < 0) ? 1 : -1)
					 * line(p[i], p[(i + 1) % p.size()]).subtendedAngle(other);
		return fabs(fabs(total) - 2.0 * PI) < EPS;
	}
	polygon halfPlaneIntersect(const line& other) const { // perform half-plane intersection with polygon. reverse line to get other half.
		vector<point> x;
		for (size_t i = 0; i < p.size(); ++i) {
			double l1 = ccw(other.a, other.b, p[i]);
			double l2 = ccw(other.a, other.b, p[(i + 1) % p.size()]);
			if (l1 > -EPS) x.push_back(p[i]);
			if (l1 * l2 < -EPS) x.push_back(other.intersectionPoint(line(p[i], p[(i + 1) % p.size()])));
		}
		return polygon(x);
	}
	double width() const { // calculate the width (smallest distance between two supporting lines) of the convex polygon
		bool first = true;
		size_t n = p.size(), tmp;
		double dist = -1.0;
		for (size_t i = 0, j = 2; ; i = (i + 1) % n) {
			line edge(p[i], p[(i + 1) % n]);
			while (edge.perpendicularDistance(p[j]) + EPS < edge.perpendicularDistance(p[(j + 1) % n])) j = (j + 1) % n;
			double cur = edge.perpendicularDistance(p[j]);
			if (!i) { // keep going until the same antipodal edge-point pair is considered
				if (first) first = false, dist = cur;
				else break;
			} else if (dist > cur) dist = cur;
			tmp = i; i = j; j = tmp;
		}
		return dist;
	}
	double diameter() const { // calculate the diameter (largest distance between two points on the polygon) of the convex polygon
		size_t n = p.size(), antipode = n;
		double dist = 0.0;
		for (size_t i = 0, j = 1; i < antipode; ++i) { // rotating calipers
			while ((p[j] - p[i]).sqnorm() + EPS < (p[(j + 1) % n] - p[i]).sqnorm()) j = (j + 1) % n;
			if (!i) antipode = j;
			dist = max(dist, (p[j] - p[i]).norm());
		}
		return dist;
	}
	double minimumDistance(const polygon& other) const { // calculate the minimum distance between two convex polygons, which can be negative.
		polygon poly[2] = {(*this), other};
		double maxseparation = 0.0;
		for (int x = 0; x < 2; ++x) {
			bool first = true;
			size_t n = poly[x].size(), no = poly[!x].size();
			double dist = 0.0;
			for (size_t i = 0, j = 0; ; i = (i + 1) % n) {
				line edge(poly[x][i], poly[x][(i + 1) % n]);
				point perp = (edge.b - edge.a).perp();
				if (i) while (perp.dot(poly[!x][j]) + EPS < perp.dot(poly[!x][(j + 1) % no])) j = (j + 1) % no;
				else {
					size_t id = 0;
					for (j = 1; j < no; ++j) if (perp.dot(poly[!x][j]) + EPS > perp.dot(poly[!x][id])) id = j;
					j = id;
				}
				double cur = (perp.dot(poly[!x][j]) - perp.dot(poly[x][i])) / perp.norm();
				if (!i) {
					if (first) first = false;
					else break;
				}
				if (cur < 0) if (dist < EPS || dist > -cur) dist = -cur;
			}
			maxseparation = x ? max(maxseparation, dist) : dist;
		}
		return maxseparation;
	}
	bool intersect(const polygon& other) const { return minimumDistance(other) < EPS; } // intersect
	double maximumDistance(const polygon& other) const { // calculate the maximum distance between two convex polygons
		bool first = true;
		size_t n = p.size(), no = other.size();
		double dist = 0.0;
		for (size_t i = 0, j = 0; ; i = (i + 1) % n) {
			line edge(p[i], p[(i + 1) % n]);
			point perp = (edge.b - edge.a).perp();
			if (i) while (perp.dot(other[j]) + EPS < perp.dot(other[(j + 1) % no])) j = (j + 1) % no;
			else {
				size_t id = 0;
				for (j = 1; j < no; ++j) if (perp.dot(other[j]) + EPS > perp.dot(other[id])) id = j;
				j = id;
			}
			double cur = max(max((p[i] - other[j]).norm(), (p[(i + 1) % n] - other[j]).norm()),
							 max((p[i] - other[(j + 1) % no]).norm(), (p[(i + 1) % n] - other[(j + 1) % no]).norm()));
			if (!i) {
				if (first) first = false;
				else break;
			}
			if (dist < cur) dist = cur;
		}
		return dist;
	}
	polygon intersection(const polygon& other) const { // finds the intersection of two convex polygons
		if (!intersect(other)) return polygon();
		polygon remainder = other;
		for (size_t i = 0; i < p.size(); ++i) {
			line scissors(p[i], p[(i + 1) % p.size()]);
			remainder = remainder.halfPlaneIntersect(scissors);
		}
		return remainder;
	}
};

// comparison by polar coordinate with reference to pivot
struct polarCompare {
	point pivot;
	polarCompare(point _pivot = point()): pivot(_pivot) {}
	bool operator()(point a, point b) {
		if (collinear(pivot, a, b)) return (a - pivot).sqnorm() < (b - pivot).sqnorm();
		return (a - pivot).angle() < (b - pivot).angle();
	}
};

// graham scan
polygon convexHull(const vector<point>& points) {
	if (points.size() <= 2) return points; // special case
	
	vector<point> p(points);
	int id = 0, n = p.size();
	for (int i = 1; i < n; ++i)
		if (p[i].y < p[id].y) id = i;
		else if (fabs(p[i].y - p[id].y) < EPS && p[i].x > p[id].x) id = i;
		
	// put pivot in front
	point tmp = p[id]; p[id] = p[0]; p[0] = tmp;
	sort(++p.begin(), p.end(), polarCompare(p[0]));
	p.push_back(p[0]); // sentinel point

	int m = 1;
	for (int i = 2; i < n + 1; ++i) {
		while (m && ccw(p[m - 1], p[m], p[i]) < EPS) --m;
		++m;
		tmp = p[m]; p[m] = p[i]; p[i] = tmp;
	}
	p.erase(p.begin() + m, p.end());
	return polygon(p);
}

// find the superposition of two vectors
void example1() {
	point a(1.0, 5.0), b(2.0, -3.0), c;
	c = a * 3.0 + b * 4.0; // note post-multiplication of doubles
	printf("%.3lf %.3lf\n", c.x, c.y);
}

// find the angle between two lines
void example2() {
	point a(0.0, 1.0), b(10.0, 3.0);
	line la(a, b), lb(5.0, 4.0, 3.0, 2.0);
	printf("%.3lf degrees\n", la.angle(lb) * 180.0 / PI);
}

// find the intersection between two line segments
void example3() {
	point a(0.0, 0.0), b(1.0, 0.0), c(1.0, 1.0), d(0.0, 1.0), res;
	line la, lb;
	
	la = line(a, b), lb = line(b, c); // contacting lines
	res = la.intersectionPoint(lb);
	if (la.contacts(lb)) printf("la contacts lb at %.3lf %.3lf\n", res.x, res.y);
	
	la = line(a, c), lb = line(b, d); // intersecting lines
	res = la.intersectionPoint(lb);
	if (la.intersects(lb)) printf("la intersects lb at %.3lf %.3lf\n", res.x, res.y);
}

// find a circle from three points on the circle = circumcircle of a triangle
void example4() {
	point a(0, 0), b(1, 0), c(1, 1);
	circle circ(a, b, c);
	printf("(%.3lf, %.3lf), %.3lf\n", circ.c.x, circ.c.y, circ.r);
}

void print(const polygon& poly) {
	printf("Polygon[");
	for (size_t i = 0; i < poly.size(); ++i) {
		printf("(%.3lf, %.3lf)", poly[i].x, poly[i].y);
		if (i < poly.size() - 1) printf(", ");
	}
	printf("]\n");
}

// perform the convex hull algorithm.
void example5() {
	vector<point> data;
	data.push_back(point(0.0, 0.0));
	data.push_back(point(5.0, 0.0));
	data.push_back(point(0.0, 10.0));
	data.push_back(point(4.0, 7.0));
	data.push_back(point(3.0, 2.0));
	data.push_back(point(6.0, 5.0));
	data.push_back(point(2.0, 8.0));
	
	// resulting hull is a counterclockwise convex polygon
	polygon hull = convexHull(data);
	printf("hull: "); print(hull);
}

// cut a square into two polygons
void example6() {
	polygon square;
	square.push_back(point(0.0, 0.0));
	square.push_back(point(1.0, 0.0));
	square.push_back(point(1.0, 1.0));
	square.push_back(point(0.0, 1.0));
	
	line scissors(0, -1, 1, 2);
	polygon left = square.halfPlaneIntersect(scissors), right = square.halfPlaneIntersect(scissors.reverse());
	printf("left: "); print(left);
	printf("right:"); print(right);
}

// compute the width (minimum distance between parallel lines touching the polygon) and diameter of a polygon (maximum distance)
void example7() {
	polygon poly;
	poly.push_back(point(1.0, 7.0));
	poly.push_back(point(0.0, 4.0));
	poly.push_back(point(0.0, 0.0));
	poly.push_back(point(1.0, 1.0));
	poly.push_back(point(2.0, 4.0));
	poly.push_back(point(2.0, 8.0));
	
	printf("poly: "); print(poly);
	
	double width = poly.width(), diameter = poly.diameter();
	printf("width %.3lf, diameter %.3lf\n", width, diameter);
}

// calculate the minimum distance / maximum distance between convex polygons
// find the intersection between convex polygons
void example8() {
	polygon polya, polyb, polyc;
	polya.push_back(point(-0.68, -4.65));
	polya.push_back(point( 2.46, -6.87));
	polya.push_back(point( 5.00, -3.61));
	polya.push_back(point( 3.04,  1.05));
	polya.push_back(point(-0.99, -0.67));
	
	// polya, translated a little so that it touches polyc
	polyb.push_back(point( 0.923, -2.383));
	polyb.push_back(point( 4.063, -4.603));
	polyb.push_back(point( 6.603, -1.343));
	polyb.push_back(point( 4.643,  3.317));
	polyb.push_back(point( 0.613,  1.597));
	
	// another polygon
	polyc.push_back(point( 2.86, 3.92));
	polyc.push_back(point(10.00, 1.50));
	polyc.push_back(point(10.00, 5.00));
	polyc.push_back(point( 5.28, 5.36));
	
	// find the intersections
	printf("polya intersect polyc %d\n", polya.intersect(polyc));
	printf("mindist %.3lf, maxdist %.3lf\n", polya.minimumDistance(polyc), polya.maximumDistance(polyc));
	printf("polyb intersect polyc %d\n", polyb.intersect(polyc));
	printf("mindist %.3lf, maxdist %.3lf\n", polyb.minimumDistance(polyc), polyb.maximumDistance(polyc));
	
	// intersection of polya and polyb = another polygon
	polygon intersection = polya.intersection(polyb);
	printf("polya: "); print(polya);
	printf("polyb: "); print(polyb);
	printf("intersect: "); print(intersection);
}

int main() {
	example1();
	example2();
	example3();
	example4();
	example5();
	example6();
	example7();
	example8();
	return 0;
}
