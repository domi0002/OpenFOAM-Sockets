/*
	C++ Socket Implementation
*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>


enum ConnectionType
{
	Server=1,
	Client=2
};

class socketFoam
{
	public:

    // Constructor depending on the connection type - server/client ?
	socketFoam(const char* hostname, int port, ConnectionType ctype);

    socketFoam(){};

	~socketFoam();

	const char* hostname;
	int port;

	template<typename E>
	ssize_t send(E* sendData, const int size);

	template<typename E>
	ssize_t recv(E* recvData, const int size);

	
	struct sockaddr_in cliServ_addr;	
	struct hostent *server;
	socklen_t clilen;
	int sockfd, newsockfd;
	ConnectionType ctype;
	void Finalize();


    void set(const char* hostname, int port, ConnectionType ctype);

	inline void Check(int ierr) { if (ierr<0) { printf("Error in communication"); exit(1) ;} }

};

template <typename E>
ssize_t socketFoam::send(E* array, const int size)
{
	if ( ctype == Client )
		return ::send(sockfd, array, size*sizeof(E),0);
	else
		return ::send(newsockfd, array, size*sizeof(E),0);
		
}

template <typename E>
ssize_t socketFoam::recv(E* array, const int size)
{
	if ( ctype == Server )
		return ::recv(newsockfd, array, size*sizeof(E),0);
	else
		return ::recv(sockfd, array, size*sizeof(E),0);
}
