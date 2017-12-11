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

#include "types.h"
#include "WavePropagation.h"
#include "scenarios/gaussian.h"
#include "writer/VtkWriter.h"
#include "tools/args.h"

#include <cstring>

int main(int argc, char** argv)
{
	// Parse command line parameters
	tools::Args args(argc, argv);

	// Scenario
	scenarios::Gaussian scenario(args.size());

	// Allocate memory across spatial domain
	Q *q = new Q[args.size()+2];

	// Initialize water height and momentum
	for (unsigned int x = 0; x < args.size()+2; x++)
	  q[x] = scenario.getQ(x);


	// Create a writer that is responsible printing out values
	//writer::ConsoleWriter writer;
	writer::VtkWriter writer("adv1d", scenario.getCellSize());

	// Helper class computing the wave propagation
	WavePropagation wavePropagation(q, args.size(), scenario.getCellSize());

	// Write initial data
	tools::Logger::logger.info("Initial data");

	// Current time of simulation
	T t = 0;

	writer.write(t, q, wavePropagation.m_maxEdgeSpeed, args.size());
	
	T maxTimeStep = args.size()/WavePropagation::ADVECTION_A * .4f;

	for (unsigned int i = 0; i < args.timeSteps(); i++) {
		// Do one time step
		tools::Logger::logger << "Computing timestep " << i
				<< " at time " << t << std::endl;

		// Update boundaries
		//		wavePropagation.setOutflowBoundaryConditions();
		wavePropagation.setPeriodicBoundaryConditions();		

		// Compute numerical flux on each edge
		maxTimeStep = wavePropagation.computeNumericalFluxes(maxTimeStep);

		// Update unknowns from net updates
		wavePropagation.updateUnknowns(maxTimeStep);

		// Update time
		t += maxTimeStep;

		// Write new values
		writer.write(t, q, wavePropagation.m_maxEdgeSpeed, args.size());
	}

	// Free allocated memory
	delete [] q;

	return 0;
}
