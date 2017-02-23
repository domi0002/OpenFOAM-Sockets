/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/


#include "socketFvPatchFieldV.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//makePatchTypeFieldTypedefs(socket);


socketFvPatchFieldV::socketFvPatchFieldV
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    socketCall(1)
{

    Info << " Constructor C" << endl;
}




socketFvPatchFieldV::socketFvPatchFieldV
(
    const socketFvPatchFieldV& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    socketCall(ptf.socketCall),
    hostPort(ptf.hostPort),
    hostAddress(ptf.hostAddress),
    communicationMode(ptf.communicationMode),
    sendSize(ptf.sendSize),
    ctype(ptf.ctype)
{
    if (&iF && mapper.hasUnmapped())
    {
        WarningIn
        (
            "socketFvPatchFieldV::socketFvPatchFieldV\n"
            "(\n"
            "    const socketFvPatchFieldV&,\n"
            "    const fvPatch&,\n"
            "    const DimensionedField<scalar, volMesh>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}



socketFvPatchFieldV::socketFvPatchFieldV
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    socketCall(1)
{
 //   evaluate();
    fvPatchField<vector>::operator=( Field<vector>("value",dict,p.size()) );
    dict.lookup("host") >> hostAddress;
    dict.lookup("port") >> hostPort;
    dict.lookup("communicationMode") >> communicationMode ;

    word cTypeStr;

    dict.lookup("connectionType") >> cTypeStr;

    if (cTypeStr=="Server")
        ctype=Server;
    else if (cTypeStr=="Client")
        ctype=Client;
    else
        FatalError << "ConnectionType must be either Server or Client " << endl;
    


    if (!Pstream::parRun() || communicationMode == "master")
    {
        Info << " host = " << hostAddress << endl;
        Info << " Port   = " << hostPort << endl ;
    }
    else
    { 
        hostPort+=Pstream::myProcNo();
        Pout << " Server = " << hostAddress << endl;
        Pout << " Port   = " << hostPort << endl ;
    }

    sendSize = this->size();

    if (communicationMode == "master" && Pstream::parRun())
    {
        sendSizeList.resize(Pstream::nProcs());
        sendSizeList[Pstream::myProcNo()]=this->size();
        Pstream::gatherList(sendSizeList);
        Pstream::scatterList(sendSizeList);

        displacements.resize(Pstream::nProcs());
        displacements[0]=0;

        for(int procID=1; procID<Pstream::nProcs(); procID++ )
        {
            displacements[procID]=displacements[procID-1]+sendSizeList[procID-1];
        }


        reduce(sendSize,sumOp<label>());
        Info <<" Send Size from OpenFOAM = " << sendSize << endl;    
        if ( Pstream::master() )
	    {
	        connector.set(hostAddress.c_str(),hostPort,ctype);
	    }

    	initializeInterpolationDataParallel();

	    computeWeights();

    }
    else // to do -- serial
    {
	//    connector.set(hostAddress.c_str(),hostPort,ctype);
    }

}

void socketFvPatchFieldV::computeWeights()
{
    vectorField delta(this->patch().delta());

    vectorField Sf(this->patch().Sf());

    scalarField magSf(this->patch().magSf());

    scalarField ndelta(this->size());
    
    // Must be changed ?
    scalarField ndeltaServer(this->size());

    commWeights.resize(this->size());

    scalarField recvBuffer;

    forAll( Sf, i )
    {
        ndelta[i]= mag((delta[i] & Sf[i])/magSf[i]);
    }

    // Server sends ndelta to client. Client interpolates it at appropriate location to ndeltaNeighbor
    List<scalarField> ndeltaGathered(Pstream::nProcs());

	// Plugin at desired locations
    ndeltaGathered[Pstream::myProcNo()]=ndelta;

    // Gather the list on the master
    Pstream::gatherList(ndeltaGathered);

    // Combine the list on master.
	scalarField  ndeltaGatheredCombined
    (
        ListListOps::combine<scalarField>
	    (
		    ndeltaGathered,
            accessOp<scalarField>()				
		)
    );

    if ( isServer() && Pstream::master())
    {
        label sendPacketSize=sendSize;
        label offset=0;
        label ierr = 0;
        while ( sendPacketSize > 0 )
        { 
       	    ierr = connector.send(ndeltaGatheredCombined.begin()+offset,sendPacketSize); 
            ierr/=sizeof(double);
            sendPacketSize-=ierr;
            offset+=ierr;
        }

    }
    else if ( isClient() && Pstream::master())
    {
        label recvPacketSize  = sendSize;
        label offset = 0;
        label ierr = 0;
        recvBuffer.resize(sendSize);

        while ( recvPacketSize > 0 )
        {
            ierr = connector.recv(recvBuffer.begin()+offset,recvPacketSize);
            ierr/=sizeof(double);
            recvPacketSize-=ierr;
            offset+=ierr;
        }

        // Client receives it and interpolates it to correct locations
        recvBuffer=mapperPtr_().interpolate(recvBuffer);
    }

    // Local parallel version
    /*
    if ( isClient() )
    {
    
        scatter(recvBuffer, ndeltaServer);

        // Note: ndeltaServer is now dimensioned by the client's dimensions

        forAll( ndeltaServer, i )
        {
            commWeights[i]=ndeltaServer[i]/(ndelta[i]+ndeltaServer[i]);
        }

        //    Pout << commWeights << endl;

    }*/


    // Root version
    if ( isClient() && Pstream::master() )
    {
        commWeights.resize(sendSize);

        forAll( commWeights, i )
        {
            commWeights[i] = recvBuffer[i]/(recvBuffer[i]+ndeltaGatheredCombined[i]);
        }
        //Info << commWeights << endl;

    }



}

void socketFvPatchFieldV::initializeInterpolationDataParallel()
{
    // Perturbation
    scalar perturb_ = 1e-5;
 	
    // Local Face Centres
    vectorField fc(this->patch().patch().faceCentres());


    //Info << this->patch().delta()<< endl;
    Info << "Weights = " << this->patch().weights()<< endl;


    // An Empty List of the gathered y coordinate
    //List<scalarField> fcGathered(Pstream::nProcs());
    List<vectorField> fcGatheredV(Pstream::nProcs());

    // Plugin at desired locations
    //fcGathered[Pstream::myProcNo()]=fcAppended;
    fcGatheredV[Pstream::myProcNo()]=fc;

    // Gather the list on the master
    //Pstream::gatherList(fcGathered);
    Pstream::gatherList(fcGatheredV);

    // Combine the list on master.
    vectorField  fcGatheredCombinedV
    (
        ListListOps::combine<vectorField>
        (
            fcGatheredV,
            accessOp<vectorField>()				
		)
    );

    scalarField fcGatheredCombined;

    if ( Pstream::master())
    {
    fcGatheredCombined.resize(3*sendSize);
    forAll( fcGatheredCombinedV, i)
    {
        fcGatheredCombined[i] = fcGatheredCombinedV[i].x();
        fcGatheredCombined[i + sendSize] = fcGatheredCombinedV[i].y();
        fcGatheredCombined[i + 2*sendSize] = fcGatheredCombinedV[i].z();
    }
    }

    scalarField samplePointsBuffer;
    vectorField samplePoints;

    if ( isServer() && Pstream::master() )
    {
        samplePointsBuffer.resize(3*sendSize);

        label sendPacketSize=sendSize*3;
        label offset=0;
        label ierr = 0;
        while ( sendPacketSize > 0 )
        { 
       	    ierr = connector.send(fcGatheredCombined.begin()+offset,sendPacketSize); 
            ierr/=sizeof(double);
            sendPacketSize-=ierr;
            offset+=ierr;
        }


        label recvPacketSize  = sendSize*3;
        offset = 0;
        while ( recvPacketSize > 0 )
        {
            ierr = connector.recv(samplePointsBuffer.begin()+offset,recvPacketSize);
            ierr/=sizeof(double);
            recvPacketSize-=ierr;
            offset+=ierr;
        }
 
    }
    else if ( isClient() && Pstream::master() )
    {  
        samplePointsBuffer.resize(3*sendSize);

        label recvPacketSize  = sendSize*3;
        label offset = 0;
        label ierr = 0;
        while ( recvPacketSize > 0 )
        {
            ierr = connector.recv(samplePointsBuffer.begin()+offset,recvPacketSize);
            ierr/=sizeof(double);
            recvPacketSize-=ierr;
            offset+=ierr;
        }

        label sendPacketSize=sendSize*3;
        offset=0;
        ierr = 0;
        while ( sendPacketSize > 0 )
        { 
       	    ierr = connector.send(fcGatheredCombined.begin()+offset,sendPacketSize); 
            ierr/=sizeof(double);
            sendPacketSize-=ierr;
            offset+=ierr;
        }
          
    }
    
    // Now exchange of faceCentres is complete on the master process. Insert into vector

    if ( Pstream::master())
    {
        samplePoints.resize(sendSize);

        forAll( samplePoints, i )
        {
            samplePoints[i].x() = samplePointsBuffer[i];
            samplePoints[i].y() = samplePointsBuffer[i + sendSize];
            samplePoints[i].z() = samplePointsBuffer[i + 2*sendSize];
        }
    }

    Pstream::scatter(samplePoints);

        // Setup interpolation .. Only valid on master() process.
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                fcGatheredCombinedV,
                perturb_
            )
        );

    // Check Interpolation
    //
    //Pout << samplePoints << endl;

   /* if ( Pstream::master() )
    {
        scalarField testField(sendSize);
        scalarField source(samplePoints.component(0));
    
        testField = mapperPtr_().interpolate(source);
        Info << "Interpolated Field =" << testField  << "Actual field = " << fcGatheredCombinedV.component(0) << endl;       
    }
    */

}

void socketFvPatchFieldV::scatter
(
    vectorField & srcv,
    vectorField & destv
)
{
    // Need step by step scatter
    for(label comp_=0;comp_<3;comp_++)
    {
        scalarField src(srcv.component(comp_));
        scalarField dest(this->size());

        MPI_Scatterv
        (
            src.begin(),
            sendSizeList.begin(),
            displacements.begin(),
            MPI_DOUBLE,
            dest.begin(),
            this->size(),
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );
        
        destv.replace(comp_,dest);
    }

}

socketFvPatchFieldV::socketFvPatchFieldV
(
    const socketFvPatchFieldV & ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    socketCall(ptf.socketCall),
    hostPort(ptf.hostPort),
    hostAddress(ptf.hostAddress),
    communicationMode(ptf.communicationMode),
    sendSize(ptf.sendSize),
    ctype(ptf.ctype)
{

    Info << " Constructor A" << endl;
}



socketFvPatchFieldV::socketFvPatchFieldV
(
    const socketFvPatchFieldV& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    socketCall(ptf.socketCall),
    hostPort(ptf.hostPort),
    hostAddress(ptf.hostAddress),
    communicationMode(ptf.communicationMode),
    sendSize(ptf.sendSize),
    ctype(ptf.ctype)
{
    Info << " Constructor B" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void socketFvPatchFieldV::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
}


void socketFvPatchFieldV::updateCoeffs()
{
   if ( updated() )
   {
	return ;
   }

    // Always send indefinitely.
    int continueSending = 1;

    int ierr=1;

    
    vectorField test;
    vectorField localBuffer(this->size());
    localBuffer = vector(0,0,0);
    

    if ( communicationMode == "master" && Pstream::parRun() )
    {
        test.resize(sendSize);
        test=vector(0,0,0);

        // Internal field to be sent over
        vectorField yc(this->patchInternalField());
	   
        // An Empty List of the gathered field
        List<vectorField> ycGathered(Pstream::nProcs());

	    // Plugin at desired locations
    	ycGathered[Pstream::myProcNo()]=yc;
	    
        // Gather the list on the master
        Pstream::gatherList(ycGathered);
	    
	    vectorField  ycGatheredCombined
    	(
            ListListOps::combine<vectorField>
	        (
		        ycGathered,
        		accessOp<vectorField>()				
		    )
        );

        if ( isServer() )
        {

            Info << "Server Vector Called " << endl;
    	    // Combine the list and send over to server
        	if ( Pstream::master() )
            {
            // Server first sends data to client
            scalarField ycSend(sendSize*3);
            for( int i = 0 ; i < sendSize ; i++)
            {
                ycSend[i]               = ycGatheredCombined[i].x();
                ycSend[i+sendSize]      = ycGatheredCombined[i].y();
                ycSend[i+2*sendSize]    = ycGatheredCombined[i].z();
            }

            label sendSizePacket=3*sendSize;
            label offset = 0;
            while ( sendSizePacket > 0 )
            {
           	    ierr = connector.send(ycSend.begin()+offset,sendSizePacket);
                ierr/=sizeof(double);
                sendSizePacket-=ierr;
                offset+=ierr; 
            }

            sendSizePacket=3*sendSize;
            offset = 0;
            while ( sendSizePacket > 0 )
            {
           	    ierr = connector.recv(ycSend.begin()+offset,sendSizePacket);
                ierr/=sizeof(double);
                sendSizePacket-=ierr;
                offset+=ierr; 
            }

            for ( int i = 0 ; i < sendSize ; i++ )
            {
                test[i].x() =   ycSend[i];
                test[i].y() =   ycSend[i+sendSize];
                test[i].z() =   ycSend[i+2*sendSize];
            }

                // Interpolate to local structure
                test=mapperPtr_().interpolate(test);


	        }
        }
        else if ( isClient() )
        {
            Info << "Client Vector Called " << endl;
            if ( Pstream::master() )
            {
                // Client first receives data
                //
                scalarField yc0(sendSize*3);
                label sendSizePacket=3*sendSize;
                label offset = 0;
                while ( sendSizePacket > 0 )
                {
           	        ierr = connector.recv(yc0.begin()+offset,sendSizePacket);
                    ierr/=sizeof(double);
                    sendSizePacket-=ierr;
                    offset+=ierr; 
                }


		        vectorField interpolatedComp(sendSize);
                forAll( interpolatedComp, i)
                {
                    interpolatedComp[i].x()=yc0[i];
                    interpolatedComp[i].y()=yc0[i+sendSize];
                    interpolatedComp[i].z()=yc0[i+2*sendSize];
                }

                //We first need to interpolate to the local structure and then surface interpolate
                interpolatedComp=mapperPtr_().interpolate(interpolatedComp);


                for ( int i = 0 ; i < sendSize ; i++ )
                {
                    yc0[i]                  = test[i].x()   =  commWeights[i]*(ycGatheredCombined[i].x()-interpolatedComp[i].x())+ interpolatedComp[i].x(); 
                    yc0[i+sendSize]         = test[i].y()   =  commWeights[i]*(ycGatheredCombined[i].y()-interpolatedComp[i].y())+ interpolatedComp[i].y(); 
                    yc0[i+2*sendSize]       = test[i].z()   =  commWeights[i]*(ycGatheredCombined[i].z()-interpolatedComp[i].z())+ interpolatedComp[i].z(); 
                }

                // Client sends data to server
                sendSizePacket=3*sendSize;
                offset = 0;
                while ( sendSizePacket > 0 )
                {
           	        ierr = connector.send(yc0.begin()+offset,sendSizePacket);
                    ierr/=sizeof(double);
                    sendSizePacket-=ierr;
                    offset+=ierr; 
                }
    
                
            }
        }
   
        scatter(test,localBuffer); 

        operator == ( localBuffer );

    }
    else // check me later
    {
        test.resize(sendSize);
        test=vector(0,0,0);

        if ( isServer() )
        {
            vectorField pIn(this->patchInternalField());
            scalarField ycSend(sendSize*3);
            for( int i = 0 ; i < sendSize ; i++)
            {
                ycSend[i]               = pIn[i].x();
                ycSend[i+sendSize]      = pIn[i].y();
                ycSend[i+2*sendSize]    = pIn[i].z();
            }

            label sendSizePacket=3*sendSize;
            label offset = 0;
            while ( sendSizePacket > 0 )
            {
           	    ierr = connector.send(ycSend.begin()+offset,sendSizePacket);
                ierr/=sizeof(double);
                Info <<" Items sent = " << ierr << endl;
                sendSizePacket-=ierr;
                offset+=ierr; 
            }

            sendSizePacket=3*sendSize;
            offset = 0;
            while ( sendSizePacket > 0 )
            {
           	    ierr = connector.recv(ycSend.begin()+offset,sendSizePacket);
                ierr/=sizeof(double);
                Info <<" Items received = " << ierr << endl;
                sendSizePacket-=ierr;
                offset+=ierr; 
            }

            for ( int i = 0 ; i < sendSize ; i++ )
            {
                test[i].x() =   ycSend[i];
                test[i].y() =   ycSend[i+sendSize];
                test[i].z() =   ycSend[i+2*sendSize];
            }

        }
        else
        {

                scalarField yc0(sendSize*3);

                vectorField pIn(this->patchInternalField());

                label sendSizePacket=3*sendSize;
                label offset = 0;
                while ( sendSizePacket > 0 )
                {
           	        ierr = connector.recv(yc0.begin()+offset,sendSizePacket);
                    ierr/=sizeof(double);
                    Info <<" Items received = " << ierr << endl;
                    sendSizePacket-=ierr;
                    offset+=ierr; 
                }

                for ( int i = 0 ; i < sendSize ; i++ )
                {
                    yc0[i]                  = test[i].x()   =   0.5*(yc0[i] + pIn[i].x());
                    yc0[i+sendSize]         = test[i].y()   =   0.5*(yc0[i+sendSize] + pIn[i].y());
                    yc0[i+2*sendSize]       = test[i].z()   =   0.5*(yc0[i+2*sendSize] + pIn[i].z());
                }

                sendSizePacket=3*sendSize;
                offset = 0;
                while ( sendSizePacket > 0 )
                {
           	        ierr = connector.send(yc0.begin()+offset,sendSizePacket);
                    ierr/=sizeof(double);
                    Info <<" Items sent = " << ierr << endl;
                    sendSizePacket-=ierr;
                    offset+=ierr; 
                }
              

        }

    operator==
    (
        test
    );


    }


    socketCall++;

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void socketFvPatchFieldV::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    
    os.writeKeyword("host") << hostAddress << token::END_STATEMENT << nl ;
    os.writeKeyword("port") <<   hostPort << token::END_STATEMENT << nl ;
    os.writeKeyword("communicationMode") << communicationMode << token::END_STATEMENT << nl;

    if ( ctype == Server )
        os.writeKeyword("connectionType") << "Server" << token::END_STATEMENT << nl;
    else if ( ctype == Client )
        os.writeKeyword("connectionType") << "Client" << token::END_STATEMENT << nl;

    this->writeEntry("value",os);

}

socketFvPatchFieldV::~socketFvPatchFieldV()
{
    //int continueSending = 0;

    if ( Pstream::parRun() && communicationMode == "master")
    {
	if ( Pstream::master() )
	{
	    //client.send(&continueSending,1);
            connector.Finalize();
	}
    }
    else
    {
        //client.send(&continueSending,1);
        //connector.Finalize();
    }
}


makePatchTypeField(fvPatchVectorField,socketFvPatchFieldV);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
