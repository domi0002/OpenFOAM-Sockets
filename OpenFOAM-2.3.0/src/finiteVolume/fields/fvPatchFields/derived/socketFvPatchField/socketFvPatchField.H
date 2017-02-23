/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::socketFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition supplies a fixed value obtained from an external 
    code using sockets

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            socketConnection;
        host            hostname; // Valid only on the client
        port            47010;
        communicationType   master; // Only root process performs the socket
        connectionType  Client;// Can be "Server" also
        value           uniform 0;
    }
    \endverbatim

Author
    Dominic Chandar
    dominic.chandar@gmail.com
    
SourceFiles
    socketFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef socketFvPatchField_H
#define socketFvPatchField_H


#include "fixedValueFvPatchFields.H"
#include "socketFoam.H"
#include "mpi.h"
#include "pointToPointPlanarInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                 Class socketFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class socketFvPatchField
:
    public fixedValueFvPatchScalarField
{
    // Private data

    label socketCall;

    label hostPort;

    word hostAddress;

    word communicationMode;

    label sendSize;

    ConnectionType ctype;

    labelList sendSizeList;

    labelList   displacements;

    scalarField commWeights;

    // 2D Planar interpolation
    autoPtr<pointToPointPlanarInterpolation> mapperPtr_;


public:

    //- Runtime type information
    TypeName("socketConnection");


    // Client or Server
    socketFoam connector;
    
    // Constructors

        //- Construct from patch and internal field
        socketFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        socketFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given socketFvPatchField
        //  onto a new patch
        socketFvPatchField
        (
            const socketFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct as copy
        socketFvPatchField
        (
            const socketFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new socketFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        socketFvPatchField
        (
            const socketFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new socketFvPatchField(*this, iF)
            );
        }




        // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );
/*
	    virtual void rmap
	    {
	 	const fvPatchScalarField&,
		const labelList&
	    }
*/

        // Check if Client or Server
        inline bool isClient(){ return ctype == Client; };

        inline bool isServer(){ return ctype == Server; };

        // Scatter utility
        void scatter(scalarField& src, scalarField& dest);

        // Initialize the interpolation between source and destination patches
        void initializeInterpolationDataParallel();
        
        // Computes weights for surface interpolation on the patch
        void computeWeights();


	    virtual void updateCoeffs();

	    ~socketFvPatchField();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "socketFvPatchField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
