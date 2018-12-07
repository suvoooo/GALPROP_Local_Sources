#ifndef Coordinate_h
#define Coordinate_h
#include "PhysicalConstants.h"
#include "pointing.h"

namespace SM {
	/** \brief Coordinates for easy access of skymaps
	 *
	 * Handles conversion between conventional l and b with the equator at b = 0 to healpix coordinates
	 */
	class Coordinate {
		private:
			double m_l, m_b; //!< The coordinates
		public:
			/**\brief Constructor that takes in galactic coordinates (l,b) in
			 * degrees.
			 *
			 * \param l the longitude in degrees, GC at 0
			 * \param b the latitude in degrees, GC at 0
			 */
			Coordinate(const double l, const double b);
			/**\brief Default constructor initializes to 0 */
			Coordinate() : m_l(0), m_b(0) {} 
			/**\brief Constructor that takes in healpix coordinate pointing */
			Coordinate(const pointing & point);

			/** \brief Return the value of l */
			const double & l() const {return m_l;}
			/** \brief Return the value of b */
			const double & b() const {return m_b;}

			/**\brief Get the angular coordinate for healpix
			 *
			 * \return a pointing object to use with healpix
			 */
			pointing healpixAng(void) const;

			/**\brief Output operator
			 *
			 * The format is (l,b)
			 */
			friend std::ostream & operator << (std::ostream &os, const Coordinate & coord);

	};
}

#endif
