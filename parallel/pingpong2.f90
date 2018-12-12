
!---------------------------------------------------------
!  Ping Pong Code
!  Based on the code on llnl.gov
!  m=0,4,16,64,256,1024,4096,16384
!  char is 1 byte.
!
!  Code is for the second part but still using nonblocking
!---------------------------------------------------------



program ping

USE mpi
integer request
integer numtasks, rank, dest, source, count, tag, ierr
integer send_stat(MPI_STATUS_SIZE)
integer stat(MPI_STATUS_SIZE)
character inmsg, outmsg

DOUBLE PRECISION start, finish
outmsg = 'x'
tag = 1



call MPI_INIT(ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

!  start = MPI_WTIME()
! task 0 sends to task 1 and waits to receive a return message

!-------------ISends and IReceives basically switch positions---------------

if (rank .eq. 0) then
   dest = numtasks-1
   source = numtasks-1

   !!!!
   call MPI_ISEND(outmsg, 0, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD,request, ierr)
   call MPI_IRECV(inmsg, 0, MPI_CHARACTER, source, tag, MPI_COMM_WORLD, stat, request, ierr)

! task 1 waits for task 0 message then returns a message
else if (rank .eq. numtasks-1) then
   dest = 0
   source = 0
   start= MPI_WTIME()
   call MPI_ISEND(outmsg, 0, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, request, ierr)
   call MPI_IRECV(inmsg, 0, MPI_CHARACTER, source, tag, MPI_COMM_WORLD, stat, request, ierr)
   finish = MPI_WTIME()
   !calls to MPI wait has been ruining the run.
   !CALL MPI_WAIT(request, send_stat, ierr)

endif

!-----------------------------------------------------------------



!finish = MPI_WTIME()


if (rank .eq. numtasks-1) then
    ! query recieve Stat variable and print message details
    write(*,*) finish-start
endif


call MPI_FINALIZE(ierr)


end 
