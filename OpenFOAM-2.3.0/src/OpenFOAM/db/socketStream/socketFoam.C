/*
	C++ Socket Implementation
*/
#include "socketFoam.H"

socketFoam::socketFoam
( 
    const char* hostname_, 
    int port_, 
    ConnectionType ctype_
)
:
hostname(hostname_),
port(port_),
ctype(ctype_)
{
	sockfd = socket(AF_INET, SOCK_STREAM, 0);
	
	if ( sockfd < 0 )
	{
		printf("Error Opening Socket\n");
		exit(1);
	}

	bzero((char *) &cliServ_addr, sizeof(cliServ_addr));

	cliServ_addr.sin_family = AF_INET;	

	if ( ctype == Server )
	{
		cliServ_addr.sin_addr.s_addr = INADDR_ANY;
		cliServ_addr.sin_port = htons(port);

        //int bufSize = 1000000;
    
        //setsockopt( sockfd, SOL_SOCKET, SO_RCVBUF, &bufSize, sizeof(int));
        //setsockopt( sockfd, SOL_SOCKET, SO_SNDBUF, &bufSize, sizeof(int));

		if ( bind( sockfd, (struct sockaddr*) &cliServ_addr, sizeof(cliServ_addr)) < 0 )
		{
			printf("Error in Binding\n");
			exit(1);
		}

		listen(sockfd, 5);

		clilen = sizeof(cliServ_addr);
		newsockfd = accept( sockfd, (struct sockaddr*) &cliServ_addr, &clilen);
		
		if ( newsockfd < 0 )
		{
			printf(" Error On Accept\n");
			exit(1);
		}

	}

	if ( ctype == Client )
	{
		server = gethostbyname(hostname);
		//bcopy((char*)server->h_addr, (char*)&cliServ_addr.sin_addr.s_addr,server->h_length);
		
		if ( inet_pton(AF_INET, hostname, &cliServ_addr.sin_addr) <= 0 )
        {
            printf("\n INET pTON error \n");
            exit(1);
        }
        else
        {
            printf("Hostname copied\n");
        }
	
		cliServ_addr.sin_port = htons(port);

		if ( connect(sockfd, (struct sockaddr *) &cliServ_addr, sizeof(cliServ_addr)) < 0 )
		{
			printf("Error in Connecting\n");
			exit(1);
		}

        
        //int bufSize = 1000000;
    
        //setsockopt( sockfd, SOL_SOCKET, SO_RCVBUF, &bufSize, sizeof(int));
    
        //setsockopt( sockfd, SOL_SOCKET, SO_SNDBUF, &bufSize, sizeof(int));

	}


}

void socketFoam::set
( 
    const char* hostname_, 
    int port_, 
    ConnectionType ctype_
)
{
    hostname = hostname_;
    port = port_;
    ctype = ctype_;

	sockfd = socket(AF_INET, SOCK_STREAM, 0);
	
	if ( sockfd < 0 )
	{
		printf("Error Opening Socket\n");
		exit(1);
	}

	bzero((char *) &cliServ_addr, sizeof(cliServ_addr));

	cliServ_addr.sin_family = AF_INET;	

	if ( ctype == Server )
	{
		cliServ_addr.sin_addr.s_addr = INADDR_ANY;
		cliServ_addr.sin_port = htons(port);

		if ( bind( sockfd, (struct sockaddr*) &cliServ_addr, sizeof(cliServ_addr)) < 0 )
		{
			printf("Error in Binding\n");
			exit(1);
		}

		listen(sockfd, 5);

		clilen = sizeof(cliServ_addr);
		newsockfd = accept( sockfd, (struct sockaddr*) &cliServ_addr, &clilen);
		
		if ( newsockfd < 0 )
		{
			printf(" Error On Accept\n");
			exit(1);
		}

	}

	if ( ctype == Client )
	{
		server = gethostbyname(hostname);
		bcopy((char*)server->h_addr, (char*)&cliServ_addr.sin_addr.s_addr,server->h_length);
		//
		
	/*	if ( inet_pton(AF_INET, hostname, &cliServ_addr.sin_addr) <= 0 )
        {
            printf("\n INET pTON error \n");
            exit(1);
        }
        else
        {
            printf("Hostname copied\n");
        }
    */
		cliServ_addr.sin_port = htons(port);

		if ( connect(sockfd, (struct sockaddr *) &cliServ_addr, sizeof(cliServ_addr)) < 0 )
		{
			printf("Error in Connecting\n");
			exit(1);
		}

	}


}
socketFoam::~socketFoam()
{


}


void socketFoam::Finalize()
{

	if (ctype==Server)
	{
		close(newsockfd);
	}

	close(sockfd);

}
