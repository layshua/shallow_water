/**
 * @file
 *  This file is part of ADV1D
 *
 *  ADV1D is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADV1D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ADV1D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Diese Datei ist Teil von ADV1D.
 *
 *  ADV1D ist Freie Software: Sie koennen es unter den Bedingungen
 *  der GNU General Public License, wie von der Free Software Foundation,
 *  Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
 *  veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
 *
 *  ADV1D wird in der Hoffnung, dass es nuetzlich sein wird, aber
 *  OHNE JEDE GEWAEHELEISTUNG, bereitgestellt; sogar ohne die implizite
 *  Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FUER EINEN BESTIMMTEN
 *  ZWECK. Siehe die GNU General Public License fuer weitere Details.
 *
 *  Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
 *  Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
 * 
 * @copyright 2017 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 * @author Leonhard Rannabauer <rannabau@in.tum.de>
 */

#include "WavePropagation.h"
#include "iostream"
#include <math.h>
T WavePropagation::computeNumericalFluxes(T dt)
{
	float maxWaveSpeed = 0.f;

	// Loop over all edges
	for (unsigned int i = 1; i < m_size+2; i++) {
		T maxEdgeSpeed;

		// Compute net updates
		LaxFriedrichsFlux(m_q[i-1], m_q[i], dt, m_cellSize,
                      m_uNetUpdatesLeft[i-1], m_uNetUpdatesRight[i-1],
                      maxEdgeSpeed );

		// Update maxWaveSpeed
		if (maxEdgeSpeed > maxWaveSpeed){
		  maxWaveSpeed = maxEdgeSpeed;
		}
	}

	// Compute CFL condition
	T maxTimeStep = m_cellSize/maxWaveSpeed * .4f;

	return maxTimeStep;
}

void WavePropagation::updateUnknowns(T dt)
{
  // Loop over all inner cells see Leveque p229 eq. (12.5)
   for (unsigned int i = 1; i < m_size+1; i++) {
     m_q[i].h -=  dt/m_cellSize * (m_uNetUpdatesRight[i-1].h + m_uNetUpdatesLeft[i-1].h);
     m_q[i].hu -=  dt/m_cellSize * (m_uNetUpdatesRight[i-1].hu + m_uNetUpdatesLeft[i-1].hu);
  }
}

void WavePropagation::setOutflowBoundaryConditions()
{
	m_q[0] = m_q[1]; 
	m_q[m_size+1] = m_q[m_size];
}

void WavePropagation::setPeriodicBoundaryConditions()
{
  m_q[m_size+1] = m_q[1];
  m_q[0] = m_q[m_size];
}



// See Leveque p. 234 eq 12.15
void WavePropagation::LaxFriedrichsFlux(Q q_l, Q q_r, T dt, T dx,
                                        Q& uNetUpdatesLeft,
                                        Q& uNetUpdatesRight,
                                        T& maxEdgeSpeed)
{

  T flux_lh  = q_l.h  * ADVECTION_A;
  T flux_lhu = q_l.hu * ADVECTION_A;
  T flux_rh  = q_r.h  * ADVECTION_A;
  T flux_rhu = q_r.hu * ADVECTION_A;
  T a = dx/dt;

  maxEdgeSpeed = fabs(ADVECTION_A);
  //  std::cout << maxEdgeSpeed << std::endl;

  uNetUpdatesRight.h  = 0.5*((flux_rh - flux_lh) - a*(q_r.h - q_l.h));
  uNetUpdatesRight.hu = 0.5*((flux_rhu - flux_lhu) - a*(q_r.hu - q_l.hu));
  uNetUpdatesLeft.h   = 0.5*((flux_rh - flux_lh) + a*(q_r.h - q_l.h));
  uNetUpdatesLeft.hu  = 0.5*((flux_rhu - flux_lhu) + a*(q_r.hu - q_l.hu));
}

