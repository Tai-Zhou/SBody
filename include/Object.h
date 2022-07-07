#ifndef SBODY_OBJECT_H
#define SBODY_OBJECT_H

#include <vector>

namespace SBody {
	namespace Object {
		class object {
		  public:
			virtual int hit(const double current[], const double last[]) = 0;
			virtual double frequency(const double current[]) = 0;
		};
		extern std::vector<object *> objectList;
		class star : public object {
		  protected:
			const int fixed;
			const double radius;
			const double r2;

		  public:
			double pos[8]; // FIXME
			star(double radius, const double position[], int fixed = 0);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double photon[]);
			double frequencyTau(const double photon[]);
		};
		class disk : public object {
		  protected:
			const double innerRadius;
			const double outerRadius;

		  public:
			disk(double innerRadius, double outerRadius);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double ph[]);
		};
		class thickDisk : public disk {
		  protected:
			const double halfAngle;

		  public:
			thickDisk(double innerRadius, double outerRadius, double halfAngle);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double ph[]);
		};
		class torus : public object {
		  protected:
			const double majorRadius;
			const double minorRadius;

		  public:
			torus(double majorRadius, double minorRadius);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double ph[]);
		};
	} // namespace Object
} // namespace SBody

#endif
