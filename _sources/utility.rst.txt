.. index::
   Utility

=======
Utility
=======

.. doxygenvariable:: sbody::absolute_accuracy
   :no-link:
.. doxygenvariable:: sbody::relative_accuracy
   :no-link:
.. doxygenvariable:: sbody::M_2PI
   :no-link:
.. doxygenvariable:: sbody::M_PI2
   :no-link:
.. doxygenvariable:: sbody::M_SQRT27
   :no-link:
.. doxygenvariable:: sbody::EPSILON
   :no-link:
.. doxygenvariable:: sbody::SIN_EPSILON
   :no-link:
.. doxygenvariable:: sbody::COS_EPSILON
   :no-link:
.. doxygenvariable:: sbody::EPSILON_CIRCLE_AREA
   :no-link:
.. doxygenvariable:: sbody::SAMPLE_NUMBER
   :no-link:
.. doxygenvariable:: sbody::ANGLE_INTERVAL
   :no-link:
.. doxygenvariable:: sbody::EPSILON_POLYGON_AREA
   :no-link:

Integrator
==========

.. doxygenclass:: sbody::Integrator
   :members:
   :no-link:

Solver
======

.. doxygenclass:: sbody::Solver
   :members:
   :no-link:


DerivativeSolver
----------------

.. doxygenclass:: sbody::DerivativeSolver
   :members:
   :no-link:


MultiSolver
===========

.. doxygenclass:: sbody::MultiSolver
   :members:
   :no-link:

MultiFunctionSolver
-------------------

.. doxygenclass:: sbody::MultiFunctionSolver
   :members:
   :no-link:

MultiDerivativeSolver
---------------------

.. doxygenclass:: sbody::MultiDerivativeSolver
   :members:
   :no-link:

GslBlock
========

.. doxygenclass:: sbody::GslBlock
   :members:
   :no-link:

Math
====

.. doxygenfunction:: sbody::Dot(const double[], const double[], size_t);
   :no-link:

.. doxygenfunction:: sbody::Dot(const double[], size_t);
   :no-link:

.. doxygenfunction:: sbody::Norm
   :no-link:

.. doxygenfunction:: sbody::Cross
   :no-link:

.. doxygenfunction:: sbody::DotCross
   :no-link:

.. doxygenfunction:: sbody::TriangleArea(double, double, double)
   :no-link:

.. doxygenfunction:: sbody::TriangleArea(const double[], const double[], const double[])
   :no-link:

.. doxygenfunction:: sbody::RotateAroundAxis
   :no-link:

.. doxygenfunction:: sbody::CartesianToSpherical(const double[], double[], size_t)
   :no-link:

.. doxygenfunction:: sbody::CartesianToSpherical(double[], size_t)
   :no-link:

.. doxygenfunction:: sbody::SphericalToCartesian(const double[], double[], size_t)
   :no-link:

.. doxygenfunction:: sbody::SphericalToCartesian(double[], size_t)
   :no-link:

.. doxygenfunction:: sbody::OppositeSign
   :no-link:

.. doxygenfunction:: sbody::MapTheta
   :no-link:

.. doxygenfunction:: sbody::ModBy2Pi
   :no-link:

.. doxygenfunction:: sbody::PhiDifference
   :no-link:

.. doxygenfunction:: sbody::LinearInterpolation(double, double, double, double, double)
   :no-link:

.. doxygenfunction:: sbody::LinearInterpolation(double, double, double, const double[], const double[], double[], size_t)
   :no-link:

.. doxygenfunction:: sbody::InterpolateSphericalPositionToCartesian
   :no-link:

.. doxygenfunction:: sbody::Flux
   :no-link:

.. doxygenfunction:: sbody::FluxDensity
   :no-link:

Elliptic Integrals
------------------

Calculate the elliptic integrals using the Legendre forms. Further information can be found in `GSL Manual <https://www.gnu.org/software/gsl/doc/html/specfunc.html#elliptic-integrals>`_, `Carlson (1988) <https://www.ams.org/mcom/1988-51-183/S0025-5718-1988-0942154-7/>`_, `Carlson (1989) <https://www.ams.org/mcom/1989-53-187/S0025-5718-1989-0969482-4/>`_, `Carlson (1991) <https://www.ams.org/mcom/1991-56-193/S0025-5718-1991-1052087-6/>`_, and `Carlson (1992) <https://www.ams.org/mcom/1992-59-199/S0025-5718-1992-1134720-4/>`_.

.. doxygenfunction:: sbody::EllipticIntegral
   :no-link:

.. doxygenfunction:: sbody::EllipticIntegral2Complex
   :no-link:

.. doxygenfunction:: sbody::EllipticIntegral4Complex
   :no-link:

.. doxygenfunction:: sbody::Carlson_RC
   :no-link:

.. doxygenfunction:: sbody::Carlson_RJ
   :no-link:
