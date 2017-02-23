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
    Foam::socketFvPatchFieldV

Group
    grpGenericBoundaryConditions

Description

    This boundary condition supplies a fixed value obtained from an external 
    code using sockets. Vector form of socketFvPatchField

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            socketConnectionV;
        host            hostname; // Valid only on the client
        port            47010;
        communicationType   master; // Only root process performs the socket
        connectionType  Client;// Can be "Server" also
        value           uniform (0 0 0);
    }
    \endverbatim

Author
    Dominic Chandar
    dominic.chandar@gmail.com


SourceFiles
    socketFvPatchFieldV.C

\*---------------------------------------------------------------------------*/

#ifndef socketFvPatchFieldV_H
#define socketFvPatchFieldV_H


#include "fixedValueFvPatchFields.H"
#include "socketFoam.H"
#include "mpi.h"
#include "pointToPointPlanarInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                 Class socketFvPatchFieldV Declaration
\*---------------------------------------------------------------------------*/
class socketFvPatchFieldV
:
    public fixedValueFvPatchVectorField
{
    // Private data

    label socketCall;

    label hostPort;

    word hostAddress;

    word communicationMode;

    label sendSize;

    ConnectionType ctype;

    labelList   sendSizeList;

    labelList   displacements;

    scalarField commWeights;

    autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

public:

    //- Runtime type information
    TypeName("socketConnectionV");


    // Client or Server
    socketFoam connector;
    
    // Constructors

        //- Construct from patch and internal field
        socketFvPatchFieldV
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        socketFvPatchFieldV
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given socketFvPatchFieldV
        //  onto a new patch
        socketFvPatchFieldV
        (
            const socketFvPatchFieldV&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct as copy
        socketFvPatchFieldV
        (
            const socketFvPatchFieldV&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new socketFvPatchFieldV(*this)
            );
        }

        //- Construct as copy setting internal field reference
        socketFvPatchFieldV
        (
            const socketFvPatchFieldV&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new socketFvPatchFieldV(*this, iF)
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
        void scatter(vectorField& src, vectorField& dest);
        
        // Initialize the interpolation between source and destination patches
        void initializeInterpolationDataParallel();

        // Computes weights for surface interpolation on the patch
        void computeWeights();


	    virtual void updateCoeffs();

	    ~socketFvPatchFieldV();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "socketFvPatchFieldV.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
