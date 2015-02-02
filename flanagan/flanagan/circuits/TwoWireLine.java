/*      Class TwoWireLine
*
*       Models a two parallel wire transmision line
*       This is a subclass of the superclass TransmissionLine
*
*       WRITTEN BY: Dr Michael Thomas Flanagan
*
*       DATE:    June 2007
*
*       DOCUMENTATION:
*       See Michael T Flanagan's Java library on-line web pages:
*       http://www.ee.ucl.ac.uk/~mflanaga/java/TransmissionLine.html
*
*       http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*       Copyright (c) June 2007    Michael Thomas Flanagan
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

public class TwoWireLine extends TransmissionLine{

    private double wireRadius = -1.0D;          // wire radius - both wires have the same radius
    private double wireSeparation = -1.0D;      // wire separation - wire centre to wire centre
    private boolean distancesSet = false;       // = true when both wire radius and separation entered

    private double relativePermittivity = 1.0D; // relative electrical permittivity of the material around the conductors

    private double relativePermeability = 1.0D; // relative magnetic permeability of the material around the conductors


    // CONSTRUCTOR
    public TwoWireLine(){
    }

    // WIRE RADII
    // Set wire radius - both wires identical
    public void setWireRadius(double radius){
        if(radius<=0.0D)throw new IllegalArgumentException("The wire radius, " + radius + ", must be greater than zero");
        if(this.wireSeparation!=-1.0D && this.wireSeparation<=2.0D*radius)throw new IllegalArgumentException("The wire separation distance, " + this.wireSeparation + ", must be greater than the sum of the two wire radii, " + 2.0D*radius);
        this.wireRadius = radius;
        if(this.wireSeparation!=-1.0)this.distancesSet = true;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // WIRE SEPARATION
    // Set wire separation - wire centre to wire centre
    public void setWireSeparation(double separation){
        if(separation<=0.0D)throw new IllegalArgumentException("The wire separation, " + separation + ", must be greater than zero");
        if(this.wireRadius!=-1.0D && separation<=2.0D*this.wireRadius)throw new IllegalArgumentException("The wire separation distance, " + separation + ", must be greater than the sum of the two wire radii, " + 2.0D*this.wireRadius);
        this.wireSeparation = separation;
        if(this.wireRadius!=-1.0)this.distancesSet = true;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PERMITTIVITY
    // Set relative electrical permittivity of the material around the conductors
    public void setRelativePermittivity(double epsilonR){
        this.relativePermittivity = epsilonR;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PERMEABILTY
    // Set relative magnetic permeability of the material around the conductors
    public void setRelativePermeability(double muR){
        this.relativePermeability = muR;
        if(this.distancesSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // CALCULATE DISTRIBUTED PARAMETERS
    private void calculateDistributedCapacitanceAndInductance(){
        super.distributedCapacitance = Impedance.parallelWiresCapacitance(1.0D, this.wireRadius, this.wireSeparation, this.relativePermittivity);
        super.distributedInductance = Impedance.parallelWiresInductance(1.0D, this.wireRadius, this.wireSeparation, this.relativePermeability);
    }
}
