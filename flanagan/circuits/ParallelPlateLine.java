/*      Class ParallelPlateLine
*
*       Models a parallel plate transmision line
*       This is a subclass of the superclass TransmissionLine
*
*       WRITTEN BY: Dr Michael Thomas Flanagan
*
*       DATE:    July 2007
*
*       DOCUMENTATION:
*       See Michael T Flanagan's Java library on-line web pages:
*       http://www.ee.ucl.ac.uk/~mflanaga/java/ParallelPlateLine.html
*       http://www.ee.ucl.ac.uk/~mflanaga/java/TransmissionLine.html
*       http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*       Copyright (c) July 2007    Michael Thomas Flanagan
*
*       PERMISSION TO COPY:
*       Permission to use, copy and modify this software and its documentation for
*       NON-COMMERCIAL purposes is granted, without fee, provided that an acknowledgement
*       to the author, Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies.
*
*       Dr Michael Thomas Flanagan makes no representations about the suitability
*       or fitness of the software for any or for a particular purpose.
*       Michael Thomas Flanagan shall not be liable for any damages suffered
*       as a result of using, modifying or distributing this software or its derivatives.
*
***************************************************************************************/

package flanagan.circuits;

import flanagan.complex.Complex;

public class ParallelPlateLine extends TransmissionLine{

    private double plateWidth = -1.0D;          // plate width
    private double plateSeparation = -1.0D;     // plate separation - inner surface to inner surface
    private boolean distancesSet = false;       // = true when both plate width and separation entered

    private double relativePermittivity = 1.0D; // relative electrical permittivity of the material between the conductors

    private double relativePermeability = 1.0D; // relative magnetic permeability of the material between the conductors


    // CONSTRUCTOR
    // Default constructor
    public ParallelPlateLine(){
        super.title = "Parallel Plate Line";
    }

    // Constructor with user suppled title
    public ParallelPlateLine(String title){
        super.title = title;
    }

    // PLATE WIDTH
    // Set plate width
    public void setPlateWidth(double width){
        if(width<=0.0D)throw new IllegalArgumentException("The plate width, " + width + ", must be greater than zero");
        this.plateWidth = width;
        if(this.plateSeparation!=-1.0)this.distancesSet = true;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PLATE SEPARATION
    // Set plate separation - inner surface to inner surface
    public void setPlateSeparation(double separation){
        if(separation<=0.0D)throw new IllegalArgumentException("The plate separation, " + separation + ", must be greater than zero");
        this.plateSeparation = separation;
        if(this.plateWidth!=-1.0)this.distancesSet = true;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PERMITTIVITY
    // Set relative electrical permittivity of the material between the conductors
    public void setRelativePermittivity(double epsilonR){
        this.relativePermittivity = epsilonR;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PERMEABILTY
    // Set relative magnetic permeability of the material between the conductors
    public void setRelativePermeability(double muR){
        this.relativePermeability = muR;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // CALCULATE DISTRIBUTED PARAMETERS
    private void calculateDistributedCapacitanceAndInductance(){
        super.distributedCapacitance = Impedance.parallelPlateCapacitance(1.0D, this.plateWidth, this.plateSeparation, this.relativePermittivity);
        super.distributedInductance = Impedance.parallelPlateInductance(1.0D, this.plateWidth, this.plateSeparation, this.relativePermeability);
    }
}
