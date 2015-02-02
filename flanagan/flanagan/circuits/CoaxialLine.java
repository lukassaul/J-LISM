/*      Class CoaxialLine
*
*       Models a coaxial transmission line
*       This is a subclass of the superclass TransmissionLine
*
*       WRITTEN BY: Dr Michael Thomas Flanagan
*
*       DATE:    July 2007
*
*       DOCUMENTATION:
*       See Michael T Flanagan's Java library on-line web pages:
*       http://www.ee.ucl.ac.uk/~mflanaga/java/CoaxialLine.html
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

public class CoaxialLine extends TransmissionLine{

    private double innerRadius = -1.0D;         // inner radius - coaxial centre to outer surface of inner conductor
    private double outerRadius = -1.0D;         // outer radius - coaxial centre to inner surface of outer conductor
    private boolean radiiSet = false;           // = true when both inner and outer radii entered

    private double relativePermittivity = 1.0D; // relative electrical permittivity of the material between the conductors

    private double relativePermeability = 1.0D; // relative magnetic permeability of the material between the conductors


    // CONSTRUCTOR
    // default constructor
    public CoaxialLine(){
        super.title = "Coaxial Line";
    }

    // user provided title
    public CoaxialLine(String title){
        super.title = title;
    }

    // RADII
    // Set inner and outer radii
    // inner: coaxial centre to outer surface of inner conductor
    // outer: coaxial centre to inner surface of outer conductor
    public void setRadii(double innerRadius, double outerRadius){
        if(innerRadius<=0.0D)throw new IllegalArgumentException("The inner radius, " + innerRadius + ", must be greater than zero");
        if(outerRadius<=0.0D)throw new IllegalArgumentException("The outer radius, " + outerRadius + ", must be greater than zero");
        if(innerRadius >= outerRadius)throw new IllegalArgumentException("The inner radius, " + innerRadius + ", must be less than the outer radius, " + outerRadius);
        this.innerRadius = innerRadius;
        this.outerRadius = outerRadius;
        this.radiiSet = true;
        this.calculateDistributedCapacitanceAndInductance();
    }

    // Set inner radius - coaxial centre to outer surface of inner conductor
    public void setInnerRadius(double radius){
        if(radius<=0.0D)throw new IllegalArgumentException("The inner radius, " + radius + ", must be greater than zero");
        if(this.outerRadius!=-1.0D && this.outerRadius<=radius)throw new IllegalArgumentException("The inner radius, " + radius + ", must be less than the outer radius, " + this.outerRadius);
        this.innerRadius = radius;
        if(this.outerRadius!=-1.0)this.radiiSet = true;
        if(this.radiiSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // Set outer radius - coaxial centre to inner surface of outer conductor
    public void setOuterRadius(double radius){
        if(radius<=0.0D)throw new IllegalArgumentException("The outer radius, " + radius + ", must be greater than zero");
        if(this.innerRadius!=-1.0D && this.innerRadius>=radius)throw new IllegalArgumentException("The outer radius, " + radius + ", must be greater than the inner radius, " + this.innerRadius);
        this.outerRadius = radius;
        if(this.innerRadius!=-1.0)this.radiiSet = true;
        if(this.radiiSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PERMITTIVITY
    // Set relative electrical permittivity of the material between the conductors
    public void setRelativePermittivity(double epsilonR){
        this.relativePermittivity = epsilonR;
        if(this.radiiSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // PERMEABILTY
    // Set relative magnetic permeability of the material between the conductors
    public void setRelativePermeability(double muR){
        this.relativePermeability = muR;
        if(this.radiiSet)this.calculateDistributedCapacitanceAndInductance();
    }

    // CALCULATE DISTRIBUTED PARAMETERS
    private void calculateDistributedCapacitanceAndInductance(){
        super.distributedCapacitance = Impedance.coaxialCapacitance(1.0D, this.innerRadius, this.outerRadius, this.relativePermittivity);
        super.distributedInductance = Impedance.coaxialInductance(1.0D, this.innerRadius, this.outerRadius, this.relativePermeability);
    }
}
