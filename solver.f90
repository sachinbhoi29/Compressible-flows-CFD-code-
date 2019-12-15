   ! Last update: 06/08/2009
   ! Higher-order version  
   ! Reference coordinates
   ! Degrees of freedom are numbered from 0 to MaxNumbDegrees
   ! Euler equations


  Module Solver

   Use basis
   Use qr
   Use Defi
   Use Mesher


  Implicit None

  
  
  ! ###############################################################################################
   Contains 
  ! ###############################################################################################


  ! ------------------------------------------------------
  Subroutine ComputeIndicatorMatrix(ir,IndicatorMatrix)
   real(8) IndicatorMatrix(MaxNumbDegrees,MaxNumbDegrees)
   integer ir,i,j,k,NumberOfPoints,StencilNum
   real(8) v1(2), v2(2),v3(2),v4(2)
   real(8) qpoints12a(2,12), qweights12a(12),qpoints12b(2,12), qweights12b(12)
   real(8) qpoints(2,25), qweights(25)
   real(8) ph, inte,check

    StencilNum = 0

       IndicatorMatrix(:,:)= 0.D0

       Select case (MeshElements(ir)%NumberOfSides)
         case(3)

          v1 = CellStencils(ir)%rCellVertexes(1,:)
          v2 = CellStencils(ir)%rCellVertexes(2,:)
          v3 = CellStencils(ir)%rCellVertexes(3,:)

          CALL Set2DIntegrationRule(qpoints12a,qweights12a,j,v1,v2,v3)
          
          qpoints(:,1:j)  = qpoints12a(:,1:j)
          qweights(1:j) = qweights12a(1:j)*TriangleArea(v1,v2,v3)
          NumberOfPoints = j

        case(4) 

          v1 = CellStencils(ir)%rCellVertexes(1,:)
          v2 = CellStencils(ir)%rCellVertexes(2,:)
          v3 = CellStencils(ir)%rCellVertexes(3,:)
          v4 = CellStencils(ir)%rCellVertexes(4,:)

          CALL Set2DIntegrationRule(qpoints12a,qweights12a,i,v1,v2,v3)
          CALL Set2DIntegrationRule(qpoints12b,qweights12b,j,v3,v4,v1)

          qpoints = 0.d0
          qpoints(:,1:i)        = qpoints12a
          qweights(1:i)         = qweights12a*TriangleArea(v1,v2,v3)
          do k=1,j
            qpoints(:,i+k+1)  = qpoints12b(:,k)
            qweights(i+k+1)   = qweights12b(k)*TriangleArea(v3,v4,v1)
          enddo


          NumberOfPoints = i+j+1

       End select ! number of sides

       ! check
!       check = 0.
!       Do i=1,NumberOfPoints
!            check = check + qweights(i)
!       Enddo 
!       check = check/LocalElementArea(StencilNum,ir,0)

        Do i=1,MaxNumbDegrees
        Do j=1,MaxNumbDegrees

       ! code for second order only  
        Inte = 0.
         Do k=1,NumberOfPoints
!           ! calculate derivatives of the reconstruction polynomial
           select case(spatialorder)
            case(2)
             ph = BasisFunctionsX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsY(qpoints(1,k),qpoints(2,k),j) 
            case(3)
             ph = BasisFunctionsX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYY(qpoints(1,k),qpoints(2,k),j)
            case(4)
             ph = BasisFunctionsX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),j)

            case(5)

             ph = BasisFunctionsX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsY(qpoints(1,k),qpoints(2,k),j) + &
    !
                  BasisFunctionsXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYYY(qpoints(1,k),qpoints(2,k),j)

            case(6)

             ph = BasisFunctionsX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsY(qpoints(1,k),qpoints(2,k),j) + &
    !
                  BasisFunctionsXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYYYY(qpoints(1,k),qpoints(2,k),j)


            case(7)

             ph = BasisFunctionsX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsY(qpoints(1,k),qpoints(2,k),j) + &
    !
                  BasisFunctionsXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYYYY(qpoints(1,k),qpoints(2,k),j) + &
!
                  BasisFunctionsXXXXXX(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXXX(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXXXY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXXY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXXYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXXYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXXYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXXYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXXYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXXYYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsXYYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsXYYYYY(qpoints(1,k),qpoints(2,k),j) + &
                  BasisFunctionsYYYYYY(qpoints(1,k),qpoints(2,k),i)*BasisFunctionsYYYYYY(qpoints(1,k),qpoints(2,k),j)


             case default 
               print*,' weno is not defined'
               read*
               stop
            end select
             Inte = Inte + ph*qweights(k)
          Enddo
            
         IndicatorMatrix(i,j) = IndicatorMatrix(i,j) + Inte
        Enddo
        Enddo

 End subroutine



! ################   Block 3: the rest  ################################################

  Subroutine  InitSolution
   Integer I,status, file_id


	! first index - RK stage
    Allocate(CSV(4,4,NElementMax), PV(4,4,NElementMax), Flux(4,4,NElementMax))
    Allocate(ExtValues(GauPointsNum,4,4,NElementMax)) ! ,ExtVertexValues(5,4,NElementMax)
    Allocate(ExtrapolatedDegreesOfFreedom(GauPointsNum,4,NelementMax,0:MaxNumbDegrees))


    Select case(SchemeType)
     case(1,2) ! linear + tvd of barth
     Allocate(DegreesOfFreedom(4,0:0,NelementMax,0:MaxNumbDegrees))
     case(3,4) ! WENO  + tvd of Tilliaeva
     Allocate(DegreesOfFreedom(4,0:4,NelementMax,0:MaxNumbDegrees))
     case default
       print*,'wrong scheme type! stop the code'
       read*
       stop
    End select

    CsV = 0.
    Pv  = 0.
    Flux = 0.
    DegreesOfFreedom = 0.
    
    ! try to read initial data
    

    file_id = 555
  
    Open(file_id,file='restart.bin',IOSTAT=STATUS,STATUS='OLD')
     IF (STATUS .EQ. 0 ) THEN    
      print*,' read the data from the restart file'
      do i=1,NElementMax
       read(file_id,*)csv(1,:,i)
      enddo
      read(file_id,*) t_
      read(file_id,*) it
      write(*,'(a20,e11.4,a24,i)')'  start from time',t_,' and time step number ',it
    
     ELSE 
    
     print*,' start from zero time'
     t_ = 0
     it = 0
     print*,' assign initial data'

    Do i=1,NElementMax 
       CSV(1,:,i) = AverageExaSol(0.D0,i)
    Enddo
    
    
!    Do i=1,NElementMax 
!     select case(IniCondType)
!     case(3)
!       if (meshelements(i)%center(1) .le. 0.1 .or. meshelements(i)%center(1) .ge. 0.2) then
!       CSV(1,:,i) = InitialData(meshelements(i)%center(1),meshelements(i)%center(2))
!      else
!       CSV(1,:,i) = AverageExaSol(0.D0,i)
!      endif      
!    case default
!!     print*,' i=',i
! !     PV(1,:,i) = InitialV(MeshElements(i)%Center(1),MeshElements(i)%Center(2))
!      CSV(1,:,i) = AverageExaSol(0.D0,i)
!    end select
!      PV(1,:,i) = CSVtoPV(CSV(1,:,i))
!    Enddo
    
    Endif ! reading from restart if - else statement
    
    ! calculate the primitive variables       
    Do i=1,NElementMax 
      PV(1,:,i) = CSVtoPV(CSV(1,:,i))
    Enddo
    
 
  End Subroutine


  ! %%%%%%%%%%%% read in the input data file & initialize everything
  Subroutine Initall
    Integer STATUS,i
    INTEGER, PARAMETER :: IOError = 111

    G1 = (gamma - 1.0)/(2.0*gamma)
    G2 = (gamma + 1.0)/(2.0*gamma)
    G3 = 2.0*gamma/(gamma - 1.0)
    G4 = 2.0/(gamma - 1.0)
    G5 = 2.0/(gamma + 1.0)
    G6 = (gamma - 1.0)/(gamma + 1.0)
    G7 = (gamma - 1.0)/2.0
    G8 = gamma - 1.0

    open(1,file='solver.ini')
     Read(1,*) IniCondType
     Read(1,*) SpatialOrder
     Read(1,*) SchemeType
     Read(1,*) limiter
     Read(1,*) FluxType
     Read(1,*) CFL
     Read(1,*) Time
     Read(1,*) StencilSizeMultiplayer
     Read(1,*) DataOutPutFreq
     Read(1,*) MovieOutPutFreq
    close(1)


   write(*,'(2x,a,2x,f8.3)') '    StencilSizeMultiplayer =',    StencilSizeMultiplayer
   write(*,*)                '    Spatial order          =',    SpatialOrder
   select case (schemetype)   
    case(1)
      write(*,*)'    Linear scheme'
    case(2)
      write(*,*)'    MUSCL'
    case(3)
      write(*,*)'    WENO'
   end select
   select case (LIMITER)   
            CASE(1)
                write(*,*)'    BARTH & JESPERSEN'
            CASE(2)
                write(*,*)'    VENKATAKRISHNAN'
   END SELECT
   write(*,*)                '    Spatial order          =',    SpatialOrder
   select case (fluxtype)   
           case(hllcfluxtype)
   write(*,*)                '    HLLC flux is used'
           case(rusanovtype)
   write(*,*)                '    Rusanov flux is used'
   end select            
   

   PolynomialOrder = SpatialOrder-1
   TheoNumDeg = (PolynomialOrder+1)*(PolynomialOrder+2)/2 - 1
   MaxNumbDegrees = TheoNumDeg


   print*
   Select case(SchemeType)
       case(1)
          print*,'Linear scheme with central stencil'
       case(2)
          print*,'MUSCL'
       case(3)
          print*,'WENO scheme'
       case(4)
          print*,'TVD scheme with Tilliaeva limiter'
   End select
   print*,' Spatial order : ',SpatialOrder
   print*,' Theoretical number of degrees of freedom = :',TheoNumDeg
   print*

  
  ! select the gaussian quadrature

 !   Xs = (xa+xb)/2  + GauPoints(j)*(xb-xa)/2

    Select case(SpatialOrder)
     case(1,2)
      GauPointsNum = 1
      GauWei(1) = 1.D0
      GauPoints(1) = 0.
      
     case(3,4)
     GauPointsNum = 2
     GauWei(1:2) = 0.5D0
     GauPoints(1) = -1./dsqrt(3.)
     GauPoints(2) = +1./dsqrt(3.)
    
    case default
     GauPointsNum = 4

     GauPoints(1)  = - Dsqrt(  (15 + dsqrt(120.))/35)
     GauPoints(2)  = - Dsqrt(  (15 - dsqrt(120.))/35)
     GauPoints(3)  = - GauPoints(2)
     GauPoints(4)  = - GauPoints(1)

     GauWei(1) = (18 - dsqrt(30.D0))/72.D0
     GauWei(2) = (18 + dsqrt(30.D0))/72.D0
     GauWei(3) = GauWei(2)
     GauWei(4) = GauWei(1)
    end select     

    ! allocate some of the arrays    

    
    
   ! initialize computational mesh
   Call InitMesh
   CALL ComputeNeibSidesID
 
   print*,' prestore reconstruction data'
   CALL     PrestoreReconstructionData
   CALL     WriteStencilInfo
   
   print*,' init solution'
   Call InitSolution
   CALL CalculateExtrapolatedDegrees

   OPEN(111,FILE = 'movie.dat')
   WRITE(111,*)'TITLE="Density" '
   WRITE(111,*)'VARIABLES="X" "Y" "D","DU","DV","E","SLOPE"'
  
   Call OutputTecplot 
  501 CONTINUE
  
  End subroutine

 ! ----------------------------------------------
 Subroutine OutPutTecplot
    integer i,k,j,counter,com,rkstage
	real(8) x,y,value(5)
  
    rkstage = 1
	write(111,*)' ZONE N=',NVerMax,' E=   ',NElementMax,'  F=FEPOINT  ET=QUADRILATERAL'
	WRITE(111,*) ', SOLUTIONTIME=',T_
	Do i=1,NVerMax
	 value = 0.0D0
	 counter = 0
	 x = Vertexes(1,i)
	 y = Vertexes(2,i)

     counter = VertexNeib(i)%iNumberOfNeib
     do k=1,counter
           j = VertexNeib(i)%iNeibIds(k)           
           value(1:4) = value(1:4) + CSV(1,1:4,j)
           VALUE(5)=VALUE(5)+MeshElements(j)%limiter
      enddo      
     value = value / counter
	 write(111,'(6(2x,e11.4))'), Vertexes(1,i), Vertexes(2,i),value(1) ,value(2),value(3),VALUE(4),value(5)
	Enddo         


    Do i=1,NElementMax
	if (MeshElements(i)%NumberOfSides .eq. 3) then
	   write(111,*) MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2), MeshElements(i)%VertexId(3), MeshElements(i)%VertexId(1)
	else
	   write(111,*) MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2), MeshElements(i)%VertexId(3), MeshElements(i)%VertexId(4)
	endif
    Enddo	       

   close(1)

 End subroutine



 ! ####################################################################
  !  here use the quadrature with  q = 6  (12 ponts)
   function IntegrateExaSolOverTriangleOld(t,v1,v2,v3)
    Real(8) IntegrateExaSolOverTriangleOld(4),t,v1(2),v2(2),v3(2),S(4)
    Real(8) qpoints(2,12), qweights(12)
    Integer k,NumberOfPoints
 
     CALL Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,v1,v2,v3)
     ! now compute the integral
     S = 0.
     Do k=1,NumberOfPoints
!          print*,'k=',k
	  S = S + InitialData(qpoints(1,k)-t,qpoints(2,k)-t)*qweights(k)
     Enddo
    ! do not forget to multiplay by the area !!
     IntegrateExaSolOverTriangleOld = S*TriangleArea(v1,v2,v3)
!     print*,' IntegrateExaSolOverTriangle =',IntegrateExaSolOverTriangle
  End function

! ####################################################################
  !  here split the triangle into four smalle ones
   function IntegrateExaSolOverTriangleLocal(t,v1,v2,v3)
    Real(8) IntegrateExaSolOverTriangleLocal(4),t,v1(2),v2(2),v3(2)
    Real(8) S1(4),S2(4),S3(4),S4(4)
    Real(8) qpoints(2,12), qweights(12)
    Integer k,NumberOfPoints
 
     ! 1sr    triangle: v1,(v1+v2)/2,(v1+v3)/2
     ! 2nd    triangle: (v1+v2)/2,(v2+v3)/2,(v1+v3)/2
     ! 3rd    triangle: (v1+v2)/2,(v2+v3)/2,v2
     ! 4th    triangle: (v1+v3)/2,(v2+v3)/2,v3     !  
  
     CALL Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,v1,(v1+v2)/2,(v1+v3)/2)
     S1 = 0.
     Do k=1,NumberOfPoints
!          print*,'k=',k
	  S1 = S1 + InitialData(qpoints(1,k)-t,qpoints(2,k)-t)*qweights(k)
     Enddo
    ! do not forget to multiplay by the area !!
     S1=  S1*TriangleArea(v1,(v1+v2)/2,(v1+v3)/2)

     CALL Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,(v1+v2)/2,(v2+v3)/2,(v1+v3)/2)
     S2 = 0.
     Do k=1,NumberOfPoints
!          print*,'k=',k
	  S2 = S2 + InitialData(qpoints(1,k)-t,qpoints(2,k)-t)*qweights(k)
     Enddo
    ! do not forget to multiplay by the area !!
     S2=  S2*TriangleArea((v1+v2)/2,(v2+v3)/2,(v1+v3)/2)

     CALL Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,(v1+v2)/2,v2,(v2+v3)/2)
     S3 = 0.
     Do k=1,NumberOfPoints
!          print*,'k=',k
	  S3 = S3 + InitialData(qpoints(1,k)-t,qpoints(2,k)-t)*qweights(k)
     Enddo
    ! do not forget to multiplay by the area !!
     S3=  S3*TriangleArea((v1+v2)/2,v2,(v2+v3)/2)

     CALL Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,(v1+v3)/2,(v2+v3)/2,v3)
     S4 = 0.
     Do k=1,NumberOfPoints
!          print*,'k=',k
	  S4 = S4 + InitialData(qpoints(1,k)-t,qpoints(2,k)-t)*qweights(k)
     Enddo
    ! do not forget to multiplay by the area !!
     S4=  S4*TriangleArea((v1+v3)/2,(v2+v3)/2,v3)

     IntegrateExaSolOverTriangleLocal  = S1+S2+S3+S4
  End function

   function IntegrateExaSolOverTriangle(t,v1,v2,v3)
    Real(8) IntegrateExaSolOverTriangle(4),t,v1(2),v2(2),v3(2)
    Real(8) S1(4),S2(4),S3(4),S4(4)
    Real(8) qpoints(2,12), qweights(12)
    Integer k,NumberOfPoints
 
!      v3
!      !   - 
!      !        -
!      !        (v2+v3)/2
!      !   -       -     -
!      !           -        -
!     (v1+v3)/2    -            v2
!      !           -        -
!      !    -      -     -
!      !      (v1+v2)/2
!      !      -
!      !  -
!      v1     
! 
     ! 1sr    triangle: v1,(v1+v2)/2,(v1+v3)/2
     ! 2nd    triangle: (v1+v2)/2,(v2+v3)/2,(v1+v3)/2
     ! 3rd    triangle: (v1+v2)/2,(v2+v3)/2,v2
     ! 4th    triangle: (v1+v3)/2,(v2+v3)/2,v3     !  
     S1 = IntegrateExaSolOverTriangleLocal(t,v1,(v1+v2)/2,(v1+v3)/2)
     S2 = IntegrateExaSolOverTriangleLocal(t,(v1+v2)/2,(v2+v3)/2,(v1+v3)/2)
     S3 = IntegrateExaSolOverTriangleLocal(t,(v1+v2)/2,v2,(v2+v3)/2)
     S4 = IntegrateExaSolOverTriangleLocal(t,(v1+v3)/2,(v2+v3)/2,v3)

     IntegrateExaSolOverTriangle  = S1+S2+S3+S4
!     print*,S1+S2+S3+S4     read*
  End function



  !  here use the quadrature with  q = 6  (12 ponts)
    function AverageExaSol(t,i)
    Real(8) AverageExaSol(4),t,v1(2),v2(2),v3(2),v4(2)
    Integer i

	! here i is the number of triangle

	Select case(MeshElements(i)%NumberOfSides)
	 case(3)
             v1 = Vertexes(:,MeshElements(i)%VertexID(1))
	     v2 = Vertexes(:,MeshElements(i)%VertexID(2))
	     v3 = Vertexes(:,MeshElements(i)%VertexID(3))
 	     AverageExaSol = IntegrateExaSolOverTriangle(t,v1,v2,v3)/MeshElements(i)%Area
     case(4)
         v1 = Vertexes(:,MeshElements(i)%VertexID(1))
	     v2 = Vertexes(:,MeshElements(i)%VertexID(2))
	     v3 = Vertexes(:,MeshElements(i)%VertexID(3))
	     v4 = Vertexes(:,MeshElements(i)%VertexID(4))

	  AverageExaSol =(IntegrateExaSolOverTriangle(t,v1,v2,v3) & 
                          +IntegrateExaSolOverTriangle(t,v3,v4,v1))/MeshElements(i)%Area
    End select
  End function

  !  -------- compute  and store side averaged basis functions value
  Subroutine  CalculateExtrapolatedDegrees
   Integer i,k,L,StencilNum,j
   Real(8) xs(2),XA(2),XB(2),XC(2)
   
   StencilNum = 0

   Do i=1,NElementMax

    Do L=1,MeshElements(i)%NumberOfSides       
          select case(MeshElements(i)%NumberOfSides)
             case(3)
               Xa =  CellStencils(i)%rCellVertexes(VertexToSide3(L,2),:)
               Xb =  CellStencils(i)%rCellVertexes(VertexToSide3(L,3),:)
             case(4)
               Xa =  CellStencils(i)%rCellVertexes(VertexToSide4(L,2),:)
               Xb =  CellStencils(i)%rCellVertexes(VertexToSide4(L,3),:)
          end select

     
         !  Gaussian rule
         do j=1,GauPointsNum
           ExtrapolatedDegreesOfFreedom(j,L,i,:) = 0.
           ExtrapolatedDegreesOfFreedom(j,L,i,0) = 1.D0

          do k=1,MaxNumbDegrees
             Xs = (xa+xb)/2  + GauPoints(j)*(xb-xa)/2
             ExtrapolatedDegreesOfFreedom(j,L,i,k) = ExtrapolatedDegreesOfFreedom(j,L,i,k) &
                 +  BasisFunctions(Xs(1),Xs(2),k,i) 
         enddo ! k
        enddo ! j 	
       Enddo  ! loop over sides
    Enddo
 
  End subroutine



  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Subroutine ComputeTimeStep(s)
   integer i,l,Rkstage
   Real(8) s, normal(2),s1,Vel,p,rad,a,b,c,diam
   Real(8) Dt_i(1:NElementMax)
   Real(8) u,v,temp,nx,ny

  
   s = 10
   Rkstage = 1
   Dt_i = 0.D0

   do i=1,NElementMax
     u = PV(Rkstage,2,i)
     v = PV(Rkstage,3,i)
     temp = PV(Rkstage,4,i)/pv(Rkstage,1,i)
     Vel = dsqrt(u**2+v**2) + dsqrt(gamma*Temp)      

     Select case (MeshElements(i)%NumberOfSides)
      case(3)
       ! for triangles use the diameter of the inner circle
       a =  MeshElements(i)%Edges(1)
       b =  MeshElements(i)%Edges(2)
       c =  MeshElements(i)%Edges(3)
       p = (a+b+c)/2
       ! radius  
!       rad = dsqrt( (p-a)*(p-b)*(p-c)/2 )
       rad = MeshElements(i)%Area/p
       MeshElements(i)%rad = RAD
       diam = 2*rad
       Dt_i(i) = diam/vel

      case(4)

        ! loop over edges
        ! use the area for time step selection
        Dt_i = 10000.D0; rad=100000000
        Do L=1,MeshElements(i)%NumberOfSides
     !     nx = MeshElements(i)%Normals(1,L)
     !     ny = MeshElements(i)%Normals(2,L)
     !     vel = abs(u*nx+v*ny)+ dsqrt(gamma*Temp)
            rad=min(rad,MeshElements(i)%Edges(L))
          Dt_i(i) = min(Dt_i(i),MeshElements(i)%Edges(L)/Vel)
        Enddo
         MeshElements(i)%rad = RAD
     !   Dt_i(i) = min(Dt_i(i),sqrt(MeshElements(i)%area)/Vel)

      End select
  enddo

  ! now select the minimun
   Do i=1,NElementMax
     s = min(s,DT_i(i))
   Enddo
   s = s*CFL

   
    if (it .le. 10 .and. IniCondType>1) s = s/10

  End subroutine

   ! Transformation from  primitive to conservative variables
   Function PVtoCSV(Pv)
    Real(8) PvToCsv(4),d,u,v,p
    REAL(8) PV(4)

    d = pv(1)
    u = pv(2)
    v = pv(3)
    p = pv(4)
    pvtocsv(1) = d
    pvtocsv(2) = d*u
    pvtocsv(3) = d*v	
    pvtocsv(4) =  p/(gamma-1) + d*(u**2 + v**2)/2
   End function


   ! Transformation from  conservative to primitive  variables
   Function CSVtoPV(CS)
    Real(8) CsvToPv(4),d,u,v,p,E
    REAL(8) CS(4)

    d = cs(1)
    u = cs(2)/cs(1)
    v = cs(3)/cs(1)
    E = cs(4)
	P = (gamma-1)*(E - d*u**2/2 - d*v**2/2)
	CSVtoPV(1) = d
	CSVtoPV(2) = u
	CSVtoPV(3) = v
	CSVtoPV(4) = p
  End function
  
 ! %%%%%%%%%%%%%%%%%%%  driver subroutine for reconstruction %%%%%%%%%%%%%%%%%%%
  Subroutine  ComputeReconstruction(Rkstage)
   integer rkstage
   real(8) s0,s1

   ! calculate degrees of freedoms
   call cpu_time(s0)
   CALL ComputeDegreesOfFreedom(Rkstage)
   call cpu_time(s1)
   LSQTime = LSQTime  + s1-s0


   call cpu_time(s0)
   select case(schemetype)
    case(1,2)
     Call ComputeReconstructionBarth(Rkstage)
    case(3)
!     print*,' no weno. press any key to stop the code'     read*     stop
    Call ComputeReconstructionWENO(Rkstage)
!    case(4)
!     Call ComputeReconstructionTVD(Rkstage)
   end select
   call cpu_time(s1)
   ExtrapTime = ExtrapTime  + s1-s0

  End subroutine
   

  ! -----------------------------------------
  SUBROUTINE ROTAFW(P, Q, COSE, SENO)
    REAL(8) COSE, P, PX, Q, QX, SENO 
       PX =  COSE*P + SENO*Q
       QX = -SENO*P + COSE*Q
       P  =  PX
       Q  =  QX
  END SUBROUTINE
 
!----------------------------------------------------------------------*
  
  SUBROUTINE ROTABK(P, Q, COSE, SENO)
    REAL(8)   COSE, P, PX, Q, QX, SENO 
      PX = COSE*P - SENO*Q
      QX = SENO*P + COSE*Q
      P  = PX
      Q  = QX
   END SUBROUTINE
 
  function flueval(CS)
    Real(8)  flueval(4)
    REAL(8)   CS(4), FLUX(4), D, U, V, P, E
      D = CS(1)
      U = CS(2)/D
      V = CS(3)/D
      E = CS(4)
      P = (gamma-1)*(E - 0.5*D*U*U- 0.5*D*V*V)

      FLUeval(1) = D*U
      FLUeval(2) = D*U*U + P
      FLUeval(3) = D*U*V
      FLUeval(4) = U*(E + P)
  end function

  function Rusanov(cdl,cdr)
    Real, dimension(4) ::  Rusanov, cdl,cdr, pvl,pvr,fdl,fdr
	Real(8) s
	 pvl = csvtopv(cdl)
	 pvr = csvtopv(cdr)
	 S = max( abs(pvl(2)) +sqrt(gamma*PVL(4)/pvl(1)),abs(pvr(2)) +sqrt(gamma*PVR(4)/pvr(1)))
	 fdl = flueval(cdl)
	 fdr = flueval(cdr)
	 Rusanov = 0.5*(fdl+fdr) - 0.5*S*(cdr-cdl)        
  end function


  

      SUBROUTINE ESTIME(PVL,PVR,SL,SR,SM)
!                                                                       
!     Purpose: to compute wave speed estimates for the HLLC Riemann solver
!  
      REAL(8)  COV,  BL,BR,CUP,GEL,GER,PM,PMAX, PMIN,PPV,PQ, PTL, PTR, QMAX, QUSER, SL, SM, SR, UM
	  REAL(8) DL,UL,VL,PL,DR,UR,VR,PR,CL,CR
	  REAL(8) PVL(4), PVR(4)
      INTEGER K
  
     !  cov - covolume, set to zero here 
      cov = 0.d0

      DL = PVL(1)
      UL = PVL(2)
      VL = PVL(3)
      PL = PVL(4)
	  CL = SQRT(Gamma*PL/DL)

      DR = PVR(1)
      UR = PVR(2)
      VR = PVR(3)
      PR = PVR(4)
	  CR = SQRT(Gamma*PR/DR)

      QUSER = 2.0
 
!     Compute guess pressure from PVRS Riemann solver
 
      CUP  = 0.25*(DL + DR)*(CL + CR)
      PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
      PPV  = MAX(0.0, PPV)
      PMIN = MIN(PL,  PR)
      PMAX = MAX(PL,  PR)
      QMAX = PMAX/PMIN
  
      IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
  
!        Select PRVS Riemann solver
 
         PM = PPV
         UM = 0.5*(UL + UR) + 0.5*(PL - PR)/CUP 
  
      ELSE
 
         BL = 1.0 - COV*DL
         BR = 1.0 - COV*DR
 
         IF(PPV.LT.PMIN)THEN
 
!           Select Two-Rarefaction Riemann solver
 
            PQ  = (PL/PR)**G1
            UM  = (PQ*UL/CL/BL + UR/CR/BR + G4*(PQ - 1.0)) 
            UM  = UM/(PQ/CL/BL + 1.0/CR/BR)
            PTL = 1.0 + G7*(UL - UM)/CL/BL
            PTR = 1.0 + G7*(UM - UR)/CR/BR
            PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
  	     
         ELSE

!           Use Two-Shock Riemann solver with PVRS as estimate
 
!           Change 20/11/2003 - introduce iterations with PVRS as initial guess
 
           DO K=1,2
             GEL = SQRT((G5*BL/DL)/(G6*PL + PPV))
             GER = SQRT((G5*BR/DR)/(G6*PR + PPV))
             PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
             UM  = 0.5*(UL + UR) + 0.5*(GER*(PM - PR) - GEL*(PM - PL))	     
             IF ( ABS((PM-PPV)/PM) .LE. 1D-4) GOTO 101
                PPV = PM
     	    ENDDO
         ENDIF
      ENDIF
      
      101 continue

!     Find speeds
 
      IF(PM.LE.PL)THEN
         SL = UL - CL
      ELSE
         SL = UL - CL*SQRT(1.0 + G2*(PM/PL - 1.0))
      ENDIF
 
      SM = UM
 
      IF(PM.LE.PR)THEN
         SR = UR + CR
      ELSE
         SR = UR + CR*SQRT(1.0 + G2*(PM/PR - 1.0))
      ENDIF
 
  END subroutine



  ! %%%%%%%%%  ! Third order TVD RK method  See Jiang G.S. and Shu C.W. (1996) and references therein
  Subroutine  ThirdOrderTVD
   Integer i,L,Rkstage
   Real(8) TotalFlux(4),s0,s1

   RkStage=1

   call cpu_time(s0)
   Call ComputeReconstruction(Rkstage)
   call cpu_time(s1)
   ReconstructionTime = ReconstructionTime + s1-s0

   call cpu_time(s0)
   Call ComputeFluxes(Rkstage)
   call cpu_time(s1)
   FluxesTime = FluxesTime + s1-s0


   call cpu_time(s0)
   Do i=1,NElementMax
     TotalFlux  =0.
     ! loop over edges
     Do L=1,MeshElements(i)%NumberOfSides
	  TotalFlux = TotalFlux +  Flux(:,L,i)*MeshElements(i)%Edges(L)
     Enddo        
     CsV(2,:,i) = CsV(1,:,i) - dt*TotalFlux/MeshElements(i)%Area
     PV(2,:,i) = CSVtoPV(CSV(2,:,i))
   Enddo 

   call cpu_time(s1)
   UPdateTime = UpdateTime  + s1-s0

   RkStage = 2

   call cpu_time(s0)
   Call ComputeReconstruction(Rkstage)
   call cpu_time(s1)
   ReconstructionTime = ReconstructionTime + s1-s0

   call cpu_time(s0)
   Call ComputeFluxes(Rkstage)
   call cpu_time(s1)
   FluxesTime = FluxesTime + s1-s0

   call cpu_time(s0)

   Do i=1,NElementMax
     TotalFlux  =0.
     ! loop over edges
     Do L=1,MeshElements(i)%NumberOfSides
	  TotalFlux = TotalFlux +  Flux(:,L,i)*MeshElements(i)%Edges(L)
     Enddo
	 CsV(3,:,i) = 0.75*CsV(1,:,i) + 0.25*CsV(2,:,i) - 0.25*dt*TotalFlux/MeshElements(i)%Area
         PV(3,:,i) = CSVtoPV(CSV(3,:,i))
   Enddo 

   call cpu_time(s1)
   UPdateTime = UpdateTime  + s1-s0


 
   RkStage = 3

   call cpu_time(s0)
   Call ComputeReconstruction(Rkstage)
   call cpu_time(s1)
   ReconstructionTime = ReconstructionTime + s1-s0

   call cpu_time(s0)
   Call ComputeFluxes(Rkstage)
   call cpu_time(s1)
   FluxesTime = FluxesTime + s1-s0

   call cpu_time(s0)

   Do i=1,NElementMax
     TotalFlux  =0.
     ! loop over edges
     Do L=1,MeshElements(i)%NumberOfSides
	  TotalFlux = TotalFlux +  Flux(:,L,i)*MeshElements(i)%Edges(L)
     Enddo
	 CsV(4,:,i) = (1./3)*CsV(1,:,i) + (2./3)*CsV(3,:,i) - (2./3)*dt*TotalFlux/MeshElements(i)%Area
         PV(4,:,i) = CSVtoPV(CSV(4,:,i))
   Enddo 

   call cpu_time(s1)
   UPdateTime = UpdateTime  + s1-s0

   CSV(1,:,:) = CSV(4,:,:)
   PV(1,:,:) = PV(4,:,:)
  end subroutine



  ! %%%%%%%%%  ! Third order TVD RK method  See Jiang G.S. and Shu C.W. (1996) and references therein
  Subroutine  FirstOrder
   Integer i,L,Rkstage
   Real(8) TotalFlux(4),s0,s1

   RkStage=1

   call cpu_time(s0)
   Call ComputeReconstruction(Rkstage)
   call cpu_time(s1)
   ReconstructionTime = ReconstructionTime + s1-s0

   call cpu_time(s0)
   Call ComputeFluxes(Rkstage)
   call cpu_time(s1)
   FluxesTime = FluxesTime + s1-s0

   call cpu_time(s0)
   Do i=1,NElementMax
     TotalFlux  =0.
     ! loop over edges
     Do L=1,MeshElements(i)%NumberOfSides
	  TotalFlux = TotalFlux +  Flux(:,L,i)*MeshElements(i)%Edges(L)
     Enddo        
!     if (i .eq. 1201) then
!       write(*,'(a,2(2x,e11.4))') 'position:',meshelements(i)%center(1),meshelements(i)%center(2)     
!       write(*,'(a,4(1x,e11.4))')' Csv:       ', CSV(1,:,i)
!       write(*,'(a,4(1x,e11.4))')' Total flux:', TotalFlux/MeshElements(i)%Area      
!       read*
!     endif
     CsV(1,:,i) = CsV(1,:,i) - dt*TotalFlux/MeshElements(i)%Area
     PV(1,:,i) = CSVtoPV(CSV(1,:,i))
     if( pv(1,1,i) < 0.   ) then
      print*,' problem in cell',i
      print*,pv(1,1,i),pv(4,1,i)
      call output
      write(*,'(a,2(2x,e11.4))') 'position:',meshelements(i)%center(1),meshelements(i)%center(2)
      TotalFlux  =0.
      ! loop over edges
      Do L=1,MeshElements(i)%NumberOfSides
	   TotalFlux = TotalFlux +  Flux(:,L,i)*MeshElements(i)%Edges(L)
      Enddo        
      Do L=1,MeshElements(i)%NumberOfSides
	    write(*,'(4(1x,e11.4))') MeshElements(i)%normals(1,L),MeshElements(i)%normals(2,L)
      Enddo        
      write(*,'(a,4(1x,e11.4))')' Total flux:', TotalFlux/MeshElements(i)%Area      
      read*
    endif      
      
   Enddo 

   call cpu_time(s1)
   UPdateTime = UpdateTime  + s1-s0
  end subroutine

 ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Subroutine ComputeConvergenceError
    real(8)  s(2), CSV0(4),dif
    integer i
  
     S = 0
     Do i=1,NelementMax
       CSV0 = AverageExaSol(time,i)
       dif  =abs( CSV0(1) - CSV(1,1,i))
       S(1) = max(S(1), dif)
       S(2) = S(2) +  dif*MeshElements(i)%Area
     Enddo
 
     open(1,file='errors.txt')
     write(1,*)' Order of accuracy:',spatialorder
     write(1,*)' Otput time:',time
     write(1,'(2x,2(a,2x,e11.4))')'   L0 error:', s(1),'   L1 error:', s(2)
     write(1,'(2x,2(2x,e11.4))')s(1),s(2)     
     write(*,'(a,2x,f11.6)')'   L0 error:', s(1)
     write(*,'(a,2x,f11.6)')'   L1 error:', s(2)
     write(1,'(4x,2(2x,f11.6))') s(1), s(2)
     close(1)
  End Subroutine 



 ! ********** Compute high-order fluxes through cell sides *********************************************************
  Subroutine ComputeFluxes(Rkstage)
   Integer i,L,j,M,norm,Rkstage,intvar,NeibId
   Real(8) xb(2),xs(2),scalar
   Real(8) normal1(2), normal2(2),s,totalflux(4)
   Real(8) DL,UL,VL,PL,DR,UR,VR,PR,COSE,SeNO
   Real(8), Dimension(4)   :: Fhllc,PVL,PVR
   Real(8), Dimension(GauPointsNum,4) :: CDL,CDR
   Integer bcnum, bctype

   ! for  the moment use only central stencil
   Flux = 0.
   Do i=1,NElementMax
     ! loop over edges
     
      Do L=1,MeshElements(i)%NumberOfSides
      
     !    L is the edge number
     
         normal1 = MeshElements(i)%Normals(:,L)
         
   	     COSE = NORMAL1(1)
         SENO = NORMAL1(2)

         CDL = 0.
         CDR = 0.

         ! 'j' is used to loop over gaussian integration points
 
!         Store initial data for LEFT state (inside cell)
         Do j=1,GauPointsNum 
           CDL(j,:) = ExtValues(j,:,L,i) 
!          Rotate forward
           CALL ROTAFW(CDL(j,2), CDL(j,3), COSE, SENO)
         Enddo

         ! Store initial data for RIGHT state

         ! find the corresponding face of the neigbour
         If (MeshElements(i)%NeighborsId(L) >0 ) then
            ! neighbour is a fluid cell
            NeibID = MeshElements(i)%NeighborsId(L) 
            Do j=1,GauPointsNum
                  intvar = GauPointsNum + 1 - j ! this is for quadrature
                  CDR(j,:) = ExtValues(intvar,:,MeshElements(i)%NeighborsSidesNum(L),NeibID) 
                  CALL ROTAFW(CDR(j,2), CDR(j,3), COSE, SENO)
             Enddo
  
         Else

          ! retrieve the number of the boundary
           BcNum = MeshElements(i)%SideBCNum(L)

           ! retrieve bc type of the boundary
           BcType  =Bctable(BcNum)%BCType
           
!           print*,bcnum,bctype

           Select case(bctype)
           case(bctypereflective)
             ! solid wall boundary
             ! print*,'reflective conditions'           
             Do j=1,GauPointsNum
                  CDR(j,1)   =   CDL(j,1)
                  CDR(j,2)   = - CDL(j,2)
                  CDR(j,3:4) =   CDL(j,3:4)
!                  print*,j                  print*,cdl(j,:)                  print*,cdr(j,:)                   
              Enddo

           case(bctypeoutflow)
           ! transmissive  boundary
             Do j=1,GauPointsNum
                  CDR(j,1)   =   CDL(j,1)
                  CDR(j,2)   = + CDL(j,2)
                  CDR(j,3:4) =   CDL(j,3:4)
              Enddo

              
              case(SYMMETRY)
              IF (IniCondType.EQ.6)THEN	!shock density interaction
                            if (MeshElements(i)%CENTer(1).lt.((1.0d0/6.0d0)+((1.0d0+20.0d0*T_)/(sqrt(3.0d0)))))then
                            PVR(1)=8.0d0
                            PVR(2)=8.25*cos(pi/6.0d0)
                            PVR(3)=-8.25*sin(pi/6.0d0)
                            PVR(4)=116.5
                            else
                            PVR(1)=1.4d0
                            PVR(2)=0.0d0
                            PVR(3)=0.0d0
                            PVR(4)=1.0d0
                            end if
                            
                            
                            
              
              end if
              
               Do j=1,GauPointsNum
                  CDR(j,:) = pvtocsv(pvr)
                  CALL ROTAFW(CDR(j,2), CDR(j,3), COSE, SENO)
                Enddo
              
              
            case(bctypefreestream)
!                print*,' inflow'
!                print*,' boundary condition zone',     MeshElements(i)%BcZone(L)
!                read*
                PVR(1)  = Bctable(bcnum)%rho
                PVR(2)  = Bctable(bcnum)%u
                PVR(3)  = Bctable(bcnum)%v
                PVR(4)  = Bctable(bcnum)%P

                Do j=1,GauPointsNum
                  CDR(j,:) = pvtocsv(pvr)
                  CALL ROTAFW(CDR(j,2), CDR(j,3), COSE, SENO)
                Enddo

            case default
             print*,' wrong bctype in the flux calculations'
             print*,bcnum,bctype             
             read*
             stop

            End select

         Endif       

        101 continue

        ! finally solve the Riemann problem at the interface
        Fhllc = 0.
        Do j=1,GauPointsNum
          select case (fluxtype)
           case(rusanovtype)
            CALL RUSANOV_FLUX(CDL(j,:),CDR(j,:),FHLLC)          
           case(hllcfluxtype)
            CALL HLLC_FLUX(CDL(j,:),CDR(j,:),FHLLC)
           case default
             print*,' flux type is not defined'
             read*
             stop            
          end select            
           !     Rotate fluxes back            
          CALL ROTABK(FHLLC(2), FHLLC(3), COSE, SENO)
          Flux(:,L,i) = Flux(:,L,i) + Fhllc*GauWei(j)
        Enddo ! j

      Enddo ! L=1,...
!
!!    if (i .eq. 1801) then
!     if (abs(meshelements(i)%center(1)-1./6) .le. 0.05 .and. meshelements(i)%center(2)>2.45) then
!       print*,'i=',i
!       write(*,'(a,2(2x,e11.4))')' Position:',meshelements(i)%center(1),meshelements(i)%center(2)     
!       print*,MeshElements(i)%NeighborsId(1:MeshElements(i)%NumberOfSides) 
!       print*
!       do l=1,MeshElements(i)%NumberOfSides
!         if (MeshElements(i)%NeighborsId(L) <0)    print*,L,MeshElements(i)%SideBCNum(L),Bctable(MeshElements(i)%SideBCNum(L))%BCType
!       enddo  
!       print*           
!       write(*,'(a,4(1x,e11.4))')' Csv:       ', CSV(1,:,i)
!       do l=1,MeshElements(i)%NumberOfSides
!       write(*,'(4(1x,e11.4))')Flux(:,L,i)
!       enddo
!       TotalFlux = 0.
!       Do L=1,MeshElements(i)%NumberOfSides
!	    TotalFlux = TotalFlux +  Flux(:,L,i)*MeshElements(i)%Edges(L)
!       Enddo        
!       write(*,'(a,4(1x,e11.4))')' Total flux:       ', Totalflux
!       read*
!     endif

   Enddo  ! i = 1, NElementMax

  End subroutine




  ! compute eigenstructure from given pair of conserved vectors
  Subroutine ComputeEigenVectors(CDL,CDR,RV,LV)
   Real(8) CDL(4),CDR(4)
   ! matrix of  right eigenvectors
   Real(8) RV(1:4,1:4)
   ! matrix of  left eigenvectors
   Real(8) LV(1:4,1:4) ! ID(1:4,1:4)
   Real(8) US,VS,VtotS,AS,HS,ES,PS,DS,S1,S2

   ! define the averaged state Us,VS,VTotS,As
   DS = 0.5*(CDL(1)+CDR(1))
   US = 0.5*(CDL(2)/CDL(1)+CDR(2)/CDR(1))
   VS = 0.5*(CDL(3)/CDL(1)+CDR(3)/CDR(1))
   ES = 0.5*(CDL(4)+CDR(4))
   ! square of the velocity
   VTOTS = US**2 + VS**2
   PS = (GAMMA-1)*(ES - 0.5*DS*VTOTS)
   AS = SQRT(GAMMA*PS/DS)
   HS = 0.5*VTOTS + AS**2/G8

   RV(1,1) =1.;        RV(1,2) = 1.;                RV(1,3) =0.;  RV(1,4) = 1.  
   RV(2,1) =Us-As;     RV(2,2) = Us;                RV(2,3) =0.;  RV(2,4) = Us+As 
   RV(3,1) =Vs;        RV(3,2) = Vs;                RV(3,3) =1.;  RV(3,4) = Vs   
   RV(4,1) =Hs-Us*As;  RV(4,2) = 0.5*VTOTS; RV(4,3) =Vs ; RV(4,4) = Hs+Us*As 

   S1 = AS/(Gamma-1)
   S2 = AS**2/(Gamma-1)

   LV(1,1) =   HS + S1*(US-AS);   LV(1,2) = -(US+S1); LV(1,3) = -VS;            LV(1,4) = 1.
   LV(2,1) =-2*HS+4*S2;           LV(2,2) = 2*US;     LV(2,3) = 2*VS;           LV(2,4) = -2.
   LV(3,1) =-2*VS*S2;             LV(3,2) =   0.;     LV(3,3) = 2*S2;           LV(3,4) = 0.
   LV(4,1) = HS- S1*(US+AS);      LV(4,2) = -US+S1;   LV(4,3) = - VS;           LV(4,4) = 1.  

   LV = LV/(2*S2)

 End subroutine

! ***************** WENO  *************************************************
  Subroutine ComputeReconstructionWENO(Rkstage)

   Integer RKStage, I,L,j
   Real(8) ExtValuesI(4,4,4),PV0(4)

   ! use only central stencil

   Do i=1,NElementMax
       ! call the reconstruction procedure for cell i
      CALL ComputeReconstructionWENO_i(Rkstage,i,ExtValuesI)
       ! verify if for any side we get negative pressure or density and adjust the data correspondingly
     Do L=1,MeshElements(i)%NumberOfSides
        Do j=1,GauPointsNum
          PV0 = csvtopv(ExtValuesI(j,:,L))
          if ( abs(pv0(1) - pv(rkstage,1,i)) > 0.9*pv(rkstage,1,i) .or. &
              abs(pv0(4) - pv(rkstage,4,i)) > 0.9*pv(rkstage,4,i) ) then
            !     use Barth limiter
              CALL ComputeReconstructionBarth_i(Rkstage,i,2,ExtValuesI)
              Goto 101
           endif
         Enddo ! j
      Enddo
     
      101 Continue

       ! finally assign the reconstructed value
       Do L=1,MeshElements(i)%NumberOfSides
          Do j=1,GauPointsNum
              ExtValues(j,:,L,i) = ExtValuesI(j,:,L)
          Enddo
       Enddo
    Enddo  ! i = 1, NElementMax

  End subroutine


! ***************** WENO  *************************************************
  Subroutine ComputeReconstructionWENO_i(Rkstage,i,ExtValuesI)
   Integer RKStage
   Integer I,k,L,j,Stmax,Lmax
   Integer StencilNum,CompNum
   Real(8) xco(2),xs(2),XA(2),XB(2),XC(2),S1,S2,S3,v1(2),v2(2),v3(2),v4(2),Sum
   Real(8) DW(4,0:4,0:MaxNumbDegrees), CharDw(4,0:4,0:MaxNumbDegrees)
   Real(8) LimitedDw(4,0:MaxNumbDegrees), FinalDw(4,0:MaxNumbDegrees),DW0(1:MaxNumbDegrees)
   Real(8) SmoothnessIndicators(4,0:MaxNumberOfStencils),ForIndicator(1:MaxNumbDegrees)
   Real(8) LinearWei(0:4),NonlinearWei(4,0:4),Tilde(4,0:4), wei
   Real(8) UnLim(4),UnLim0(4), UnLim1(4),UnLim2(4),UnLim3(4)  
   Real(8) ExtValuesI(4,4,4)

   ! eigenstructure
   Real(8) CDL(4),CDR(4)
   ! matrix of  right eigenvectors
   Real(8) RV(1:4,1:4)
   ! matrix of  left eigenvectors
   Real(8) LV(1:4,1:4) ! ID(1:4,1:4)
   Real(8) Ie(4,4) 
   Real(8) normal1(2),COSE,SeNO
   Real(8) PV0(4)
   real::eps_weno,power
   Integer bcnum,bctype
    eps_weno=1e-6
    power=4
 
    Do L=1,MeshElements(i)%NumberOfSides       
   
       CDL = CSV(RKSTAGE,:,i)
       If (MeshElements(i)%NeighborsId(L) >0) then
       
         CDR = CSV(RKSTAGE,:,MeshElements(i)%NeighborsId(L))
         
       Else
          ! retrieve the number of the boundary
           BcNum = MeshElements(i)%SideBCNum(L)
           ! retrieve bc type of the boundary
           BcType  =Bctable(BcNum)%BCType
         ! solid boundary conditions
          Select case(bctype)
            case(bctypereflective)
             CDR = CDL
             CDR(2)   = - CDL(2)
            case default
           ! transmissive  boundary
             CDR   =   CDL
            End select
       Endif

       normal1 = MeshElements(i)%Normals(:,L)
       COSE = NORMAL1(1)
       SENO = NORMAL1(2)
    
       CALL ROTAFW(CDL(2), CDL(3), COSE, SENO)
       CALL ROTAFW(CDR(2), CDR(3), COSE, SENO)

       Call ComputeEigenVectors(CDL,CDR,RV,LV)

       ! initialize linear weights 
       LinearWei(0)   = 100.   !central stencil linear weight (possible choices are 10,100,1000)
       LinearWei(1:CellStencils(i)%iNumberOfStencils) = 1.
       Tilde = 0.
       do compnum=1,4
         Tilde(compnum,0:CellStencils(i)%iNumberOfStencils) = LinearWei(0:CellStencils(i)%iNumberOfStencils)
       enddo

       ! retrieve degrees of freedom for all stencils & components of the conserved vector
       DW(:,0:CellStencils(i)%iNumberOfStencils,:) = DegreesOfFreedom(:,0:CellStencils(i)%iNumberOfStencils,i,:)

       ! project all degrees of freedom on normal direction 
       CharDW = 0.D0  

       Do StencilNum=0,CellStencils(i)%iNumberOfStencils
       Do k=0,MaxNumbDegrees
        CharDw(:,StencilNum,k) = matmul(lv,Dw(:,StencilNum,k))
       Enddo
       Enddo


      ! Loop over all stencils & components of the conserved vector
      ! compute and store the smoothness indicators for each stencil & each component
       Do StencilNum=0,CellStencils(i)%iNumberOfStencils
         Do CompNum=1,4
           DW0 = CharDw(CompNum,StencilNum,1:MaxNumbDegrees)
           ForIndicator = matmul(CellStencils(i)%rIndicatorMatrix(:,:),Dw0)
           SmoothnessIndicators(CompNum,StencilNum) = dot_product(Dw0,Forindicator)
         Enddo
       Enddo

      ! compute the final extrapolated value for each side L

      ! nonlinear weights
       do compnum=1,4
        do stencilnum=0,CellStencils(i)%iNumberOfStencils
         Tilde(compnum,stencilnum) = LinearWei(stencilnum)/(eps_weno + SmoothnessIndicators(compnum,stencilnum))**power
        enddo
       enddo

       do compnum=1,4
         Sum = 0.
         Do stencilnum=0,CellStencils(i)%iNumberOfStencils
            Sum = Sum + Tilde(compnum,stencilnum)
         Enddo
         Do StencilNum=0,CellStencils(i)%iNumberOfStencils
           NonlinearWei(CompNum,StencilNum) = Tilde(CompNum,StencilNum)/Sum
         Enddo
      enddo

      ! limit the degrees of freedom:
      ! DW(4,0:4,0:MaxNumbDegrees)

       LimitedDw(:,0:MaxNumbDegrees) = 0.D0
       Do k=0,MaxNumbDegrees
         Do StencilNum=0,CellStencils(i)%iNumberOfStencils
          LimitedDw(:,k) = LimitedDw(:,k) + CharDW(:,StencilNum,k)* NonlinearWei(:,StencilNum)
         Enddo
       Enddo

       ! project back
       Do k=0,MaxNumbDegrees
          FinalDw(:,k) = Matmul(RV,LimitedDw(:,k))
       Enddo

!       FinalDW = LimitedDW

       ! finally compute WENO reconstructed value for side L
       StencilNum = 0
       select case(MeshElements(i)%NumberOfSides)
          case(3)
            Xa =  CellStencils(i)%rCellVertexes(VertexToSide3(L,2),:)
            Xb =  CellStencils(i)%rCellVertexes(VertexToSide3(L,3),:)
          case(4)
            Xa =  CellStencils(i)%rCellVertexes(VertexToSide4(L,2),:)
            Xb =  CellStencils(i)%rCellVertexes(VertexToSide4(L,3),:)
       end select

       do j=1,GauPointsNum
           Xs = (xa+xb)/2  + GauPoints(j)*(xb-xa)/2
           Unlim = FinalDw(:,0) 
           do k=1,MaxNumbDegrees
                wei = ExtrapolatedDegreesOfFreedom(j,L,i,k)
                Unlim = Unlim +  wei*FinalDw(:,k)
           enddo ! k
           ExtValuesI(j,:,L) =  Unlim
       enddo ! j

       Enddo  ! loop over sides

  End subroutine


! ***************** Linear & TVD of Barth  *************************************************
  Subroutine ComputeReconstructionBarth(Rkstage)
  
   Integer RKStage, I,L,j
   Real(8) ExtValuesI(4,4,4),PV0(4)

   ! use only central stencil

   Do i=1,NElementMax
       ! call the reconstruction procedure for cell i
       CALL ComputeReconstructionBarth_i(Rkstage,i,schemetype,ExtValuesI)
       ! finally assign the reconstructed value
       Do L=1,MeshElements(i)%NumberOfSides
          Do j=1,GauPointsNum
              ExtValues(j,:,L,i) = ExtValuesI(j,:,L)
          Enddo
       Enddo
    Enddo  ! i = 1, NElementMax

  End subroutine



! ***************** Linear & TVD of Barth  *************************************************
  Subroutine ComputeReconstructionBarth_i(Rkstage,i,schemetype,ExtValuesI)
   IMPLICIT NONE
   Integer RKStage,schemetype
   Integer I,k,L,j
   Integer StencilNum 
   Real(8) xco(2),xs(2),XA(2),XB(2),XC(2),S1,S2,S3,v1(2),v2(2),v3(2),v4(2),S
   Real(8) DW(0:MaxNumbDegrees),US,Unlim, Ulimited
   ! for the limiter
   Real(8) VerL(2,5), Vert(5),Vi(5),PV0(4), HighOrderTerms(4)
   Real(8) VA,PSI(4),UMIN,UMAX,SJ,limiter2
   Integer NeID,CompNum,NumOfNeighLim
   Real(8) ExtValuesI(4,4,4)
   real::DMIN,DPLUS,SFD,D2,KAPPA_VEN,DELTA_X,EPSI2,EPSI_ZERO,y


   ! use only central stencil

      StencilNum = 0
     ! loop over components of the conserved vector
      Do CompNum=1,4
       ! 1. retrieve degrees of freedom from central stencil
       DW = DegreesOfFreedom(CompNum,StencilNum,i,:)

       ! 2. compute the limiter psi and also count the number of neigbouring cells used in the limiting procedure
      Select case(schemetype)

      case default
       print*,' wrong scheme type in barth procedure'
       read*
       stop

      case(1) ! linear scheme
       limiter2 = 1.D0

       case(2) ! TVD scheme
       !write(400,*)LIMITER,"KOITA EDW"
       NumOfNeighLim = 0
       Do L=1,MeshElements(i)%NumberOfSides
          NeId = MeshElements(i)%NeighborsID(L)
          if (NeID <0) then
!           print*,' ups'
            Cycle
          endif
          NumOfNeighLim = NumOfNeighLim + 1
          VerL(:,NumOfNeighLim) = CellStencils(i)%Stencils(StencilNum)%rElementCentres(L,:) 
          Vi(NumOfNeighLim) = CSV(Rkstage,CompNum,NeID)
      Enddo ! L =1

      psi = 0.
      
      Va = CSV(Rkstage,CompNum,i)

      umin = Va
      umax = Va
	  
      Do L=1,NumOfNeighLim
	   umin = min(umin,Vi(L))
	   umax = max(umax,Vi(L))
      Enddo

	  ! compute the values at the vertices
      Vert = 0.
      HighOrderTerms = 0.
      Do L=1,MeshElements(i)%NumberOfSides

!         Select case(MeshElements(i)%NumberOfSides)
!             case(3)
!               Xa =  CellStencils(i)%rCellVertexes(VertexToSide3(L,2),:)
!               Xb =  CellStencils(i)%rCellVertexes(VertexToSide3(L,3),:)
!  	           Xs = (Xa+Xb)/2
!             case(4)
!               Xa =  CellStencils(i)%rCellVertexes(VertexToSide4(L,2),:)
!               Xb =  CellStencils(i)%rCellVertexes(VertexToSide4(L,3),:)
!  	           Xs = (Xa+Xb)/2
!         End select

        !   most restrictive - check vertexes but only first order for triangles
          Xs = CellStencils(i)%rCellVertexes(L,:) !  LocalElementVertexes(StencilNum,i,0,L,:)

         ! first store higher-order terms
         do k=1,MaxNumbDegrees
            HighOrderTerms(L) = HighOrderTerms(L) +  Basisfunctions(xs(1),xs(2),k,i)*dw(k)
         enddo ! k
        Vert(L) = DW(0) + HighOrderTerms(L)
      Enddo

      
	  ! compute psi
  Do L=1,MeshElements(i)%NumberOfSides
      
            IF (LIMITER.EQ.1)THEN !BARTH AND JESPERSEN
	   if ( abs(Vert(L)-Va) .le. 1e-6) then
	    psi(L) = 1.
	    
           else
	     if (Vert(L)-Va >0.) then
              psi(L) = min(1., (umax-va)/(vert(L)-va))
              
             else
              psi(L) = min(1., (umin-va)/(vert(L)-va))
              
             endif
	   endif
            END IF
	   IF (LIMITER.EQ.2)THEN !VENKATATKRISHNAN
	   D2=(Vert(L)-Va);KAPPA_Ven=10.0D0 ;DELTA_X=MeshElements(i)%rad
	   EPSI_ZERO=1E-16
	   
                IF (abs(d2).le. 1e-6)THEN
                    psi(L) = 1.0D0
                ELSE
                        IF (D2.GT.0.0D0)THEN
                              DMIN=D2
                              DPLUS=(umax-va)
                              EPSI2=(KAPPA_VEN*MeshElements(i)%rad)**3
                              psi(L)=(1/DMIN)*(((DPLUS**2+EPSI2)*DMIN+2*(DMIN**2)*DPLUS)/(DPLUS**2+2*DMIN**2+DMIN*DPLUS+EPSI2))
                
                        END IF
                        IF (D2.LT.0.0D0)THEN

			      DMIN=D2
                              DPLUS=(umin-va)
                              EPSI2=(KAPPA_VEN*MeshElements(i)%rad)**3
                              psi(L)=(1/DMIN)*(((DPLUS**2+EPSI2)*DMIN+2*(DMIN**2)*DPLUS)/(DPLUS**2+2*DMIN**2+DMIN*DPLUS+EPSI2))
                    
                          !y=(umin-va)/(vert(L)-va)
                    
                          !PSI(L)=min(1., (y**2+2*y)/(y**2+y+2))
                    
                    
                        END IF
                END IF
          END IF
          
	  End do
      
      
    IF  (CompNum.EQ.1)THEN
        MeshElements(i)%LIMITER=limiter2
        
    END IF
    
     End select ! scheme type & calculation of the limiter

      ! 3.  compute  and store  reconstructed value
        Do L=1,MeshElements(i)%NumberOfSides       
          select case(MeshElements(i)%NumberOfSides)
             case(3)
               Xa =  CellStencils(i)%rCellVertexes(VertexToSide3(L,2),:)
               Xb =  CellStencils(i)%rCellVertexes(VertexToSide3(L,3),:)
             case(4)
               Xa =  CellStencils(i)%rCellVertexes(VertexToSide4(L,2),:)
               Xb =  CellStencils(i)%rCellVertexes(VertexToSide4(L,3),:)
          end select
     
          !  store extrapolated values in gaussian points
         do j=1,GauPointsNum
           Xs = (xa+xb)/2  + GauPoints(j)*(xb-xa)/2
           Unlim = DW(0) ! CSV(Rkstage,CompNum,i)
           do k=1,MaxNumbDegrees
                Unlim = Unlim +  limiter2*ExtrapolatedDegreesOfFreedom(j,L,i,k)*dw(k)
           enddo ! k

           ExtValuesI(j,CompNum,L) =  Unlim

        enddo ! loop over gaussian integration points j 	

       Enddo  ! loop over sides

     Enddo ! loop over components

      ! verify if for any side we get negative pressure or density and adjust the data correspondingly
     Do L=1,MeshElements(i)%NumberOfSides
        Do j=1,GauPointsNum
          PV0 = csvtopv(ExtValuesI(j,:,L))
          if ( abs(pv0(1) - pv(rkstage,1,i)) > 0.9*pv(rkstage,1,i) .or. &
              abs(pv0(4) - pv(rkstage,4,i)) > 0.9*pv(rkstage,4,i) ) then
            !     use first order'
              ExtValuesI(j,:,L) = CSV(rkstage,:,i)
           endif
         Enddo ! j
      Enddo
 


  End subroutine



 !%%%%%%%%%%%%%%%%%%%%%%%%%%% Output procedure %%%%%%%%%%%%%%%5


  ! -------------------------------------------
  Subroutine Output
   integer i,k,j,counter,com,rkstage
   real(8) value(4), Xs(2)

   rkstage = 1

   
 !  print*,' calculate degrees of freedom'
!   CALL ComputeDegreesOfFreedom(Rkstage)

 !  print*,' calculate reconstruction values in vertexes'
!   CALL ComputeReconstructionVertexes

    open(1,file='rho.dat')
    write(1,*) 'TITLE = "rho""'
    write(1,*)' VARIABLES = "x" "y" "rho"'
   
    open(2,file='p.dat')
    write(2,*) 'TITLE = "p""'
    write(2,*)' VARIABLES = "x" "y" "p"'

    open(3,file='vel.dat')
    write(3,*) 'TITLE = "vel""'
    write(3,*)' VARIABLES = "x" "y" "vel"'

    open(4,file='stream.dat')
    write(4,*) 'TITLE = "stream""'
    write(4,*)' VARIABLES = "x" "y" "u" "v"'

!    open(5,file='vorticity.dat')
!    write(5,*) 'TITLE = "vorticity""'
!    write(5,*)' VARIABLES = "x" "y" "vor"'

    do k=1,4
        write(k,*)' ZONE N=',NVerMax,' E=   ',NElementMax,'  F=FEPOINT  ET=QUADRILATERAL'
         WRITE(k,*) ', SOLUTIONTIME=',T_
    enddo

!   print*,' start output'

    Do i=1,NVerMax
         value = 0.
         counter = 0
         ! vertex coordinate
         XS = Vertexes(:,i)

         counter = VertexNeib(i)%iNumberOfNeib
!        print*,counter         read*
         do k=1,counter
!           print*,VertexNeib(i)%iNeibIds(k)           
           j = VertexNeib(i)%iNeibIds(k)           
           value = value + PV(1,:,j)             
         enddo      
         value = value / counter

         write(1,'(5(2x,e11.4))'), Vertexes(1,i), Vertexes(2,i),value(1)
         write(2,'(5(2x,e11.4))'), Vertexes(1,i), Vertexes(2,i),value(4)
         write(3,'(5(2x,e11.4))'), Vertexes(1,i), Vertexes(2,i),dsqrt( value(2)**2+value(3)**2)
         write(4,'(5(2x,e11.4))'), Vertexes(1,i), Vertexes(2,i),value(2),value(3)
!         write(5,'(5(2x,e11.4))'), Vertexes(1,i), Vertexes(2,i),abs(value(5))
    Enddo

   
    Do i=1,NElementMax
       do k=1,4
        if (MeshElements(i)%NumberOfSides .eq. 3) then
           write(k,*) MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2),  MeshElements(i)%VertexId(3),  MeshElements(i)%VertexId(1)
   	    else
           write(k,*) MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2),  MeshElements(i)%VertexId(3),  MeshElements(i)%VertexId(4)
        endif
       enddo                 
    Enddo	  
     
    close(1)
    close(2)     
    close(3)
    close(4) 
!    close(5)

  End subroutine






 ! ############# Compute degrees of freedom for each stencil
  Subroutine ComputeDegreesOfFreedom(Rkstage)

   Integer Rkstage,CompNum
   Integer I,j,Stmax, StencilNum 
   Real(8), allocatable :: Subst(:)
   Real(8) DW(1:MaxNumbDegrees)


     Do i=1,NElementMax ! loop over spatial cells
      
      ! Loop over all stencils
      Do StencilNum=0,CellStencils(i)%iNumberOfStencils   !	  MeshElements(i)%NumberOfStencils

        StMax = CellStencils(i)%Stencils(StencilNum)%iSize

		allocate(Subst(Stmax))

        Do CompNum=1,4 ! loop over components of conserved vector

         DegreesOfFreedom(CompNum,StencilNum,i,:) = 0.D0

   !        calculate the right hand side
	  	 Subst = 0.
		 do j=1,StMax
		  Subst(j) = (Csv(Rkstage,Compnum,CellStencils(i)%Stencils(StencilNum)%iCellIDs(j)) - &
                     Csv(Rkstage,Compnum,i)) 
             ! V(Rkstage,CellStencils(i)%Stencils(StencilNum)%iCellIDs(j)) - V(Rkstage,i)
		 enddo

          ! calculate degrees of freedom
  	    Dw(1:MaxNumbDegrees) = matmul(CellStencils(i)%Stencils(StencilNum)%rFinalMatrix,Subst)

   !    print*,' 2. store the degrees'
        DegreesOfFreedom(CompNum,StencilNum,i,0) = CSV(rkstage,compnum,i)
        DegreesOfFreedom(CompNum,StencilNum,i,1:MaxNumbDegrees) = DW
       Enddo ! compnum

	   deallocate(Subst)

      Enddo ! stencil num

    Enddo ! i

  End subroutine
  
  

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Subroutine WriteStencilInfo
   integer i,l,file_id,k
 
   101 format(2x,a,2x,i4)
   201 format(2x,a,2(2x,e11.4))

   file_id = 555
   open(file_id,file='stencilinfo.txt')
    
   write(file_id,'(2x,a,2x,f8.3)') '    StencilSizeMultiplayer =',    StencilSizeMultiplayer
   write(file_id,*)                '    Spatial order          =',    SpatialOrder
   select case (schemetype)   
    case(1)
      write(file_id,*)'    Linear scheme'
    case(2)
      write(file_id,*)'    TVD of Barth'
    case(3)
      write(file_id,*)'    WENO'
   end select
       
   write(file_id,*)                '    Spatial order          =',    SpatialOrder
   select case (fluxtype)   
           case(hllcfluxtype)
   write(file_id,*)                '    HLLC flux is used'
           case(RUSANOVtype)
   write(file_id,*)                '    RUSANOV flux is used'
   end select            
   
   Write(file_id,*) ' Boundary condition types'

   do k=1,NBSets
      write(file_id,*) ' --------------'
      write(file_id,*) '   ',Bctable(k)%SurfaceName
      write(file_id,*) Bctable(k)%BcType
      write(file_id,*) Bctable(k)%rho
      write(file_id,*) Bctable(k)%u
      write(file_id,*) Bctable(k)%v
      write(file_id,*) Bctable(k)%P
   enddo
   write(file_id,*) ' --------------'

   write(file_id,*) ''
   WRITE(file_id,'(a,i4)')' Total number of cells ', NElementMax
   WRITE(file_id,'(a,i4)')' Total number of triangles ', NTriMax
   WRITE(file_id,'(a,i4)')' Total number of quadrilaterals ', NQuadMax
   
   
   Do i = 1, NElementMax
   write(file_id,*)' ----------------------------------- '
   write(file_id,101)' Cell  number                ',i 
   write(file_id,101)' Number of sides of  the cell:',CellStencils(i)%iNumberOfSides
   write(file_id,101)' Number of stencils          :',CellStencils(i)%iNumberOfStencils
   write(file_id,201)' Cell centre coordinates     :',MeshElements(i)%Center
   Do l=0,CellStencils(i)%iNumberOfStencils
     write(file_id,101)'    Stencil number:',L
     write(file_id,101)'    Number of cells in the stencil:',CellStencils(i)%Stencils(L)%iSize
   Enddo
   write(file_id,*)''
   Enddo
   close(file_id)

  End subroutine

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Subroutine WriteInfo
   integer i,l,file_id,k
   Real TotalMemory
   
 
   101 format(2x,a,2x,i4)
   201 format(2x,a,2(2x,e11.4))

   file_id = 555
   open(file_id,file='info.txt')
    
   write(file_id,'(2x,a,2x,f8.3)') '    StencilSizeMultiplayer =',    StencilSizeMultiplayer
   write(file_id,*)                '    Spatial order          =',    SpatialOrder
   select case (schemetype)   
    case(1)
      write(file_id,*)'    Linear scheme'
    case(2)
      write(file_id,*)'    TVD of Barth'
    case(3)
      write(file_id,*)'    WENO'
   end select    
       
   write(file_id,*)                '    Spatial order          =',    SpatialOrder
   select case (fluxtype)   
           case(hllcfluxtype)
   write(file_id,*)                '    HLLC flux is used'
           case(RUSANOVtype)
   write(file_id,*)                '    RUSANOV flux is used'
   end select            
   
   Write(file_id,*) ' Boundary condition types'

   do k=1,NBSets
      write(file_id,*) ' --------------'
      write(file_id,*) '   ',Bctable(k)%SurfaceName
      write(file_id,*) Bctable(k)%BcType
      write(file_id,*) Bctable(k)%rho
      write(file_id,*) Bctable(k)%u
      write(file_id,*) Bctable(k)%v
      write(file_id,*) Bctable(k)%P
   enddo
   write(file_id,*) ' --------------'

   write(file_id,*) ''
   WRITE(file_id,'(a,i)')' Total number of cells ', NElementMax
   WRITE(file_id,'(a,i)')' Total number of triangles ', NTriMax
   WRITE(file_id,'(a,i)')' Total number of quadrilaterals ', NQuadMax   

   TotalMemory =   size(csv) + size(pv)  +  size(flux)  + size(DegreesOfFreedom) +size(ExtrapolatedDegreesOfFreedom)  &
                  + size(ExtValues) + size(ExtVertexValues)  + size(cellstencils)    + size(meshelements)

   222 format(2x,a,1x,e11.4)
   write(file_id,'(2x,a,1x,i)') 'NElementMax = ',NElementMax
   write(file_id,222) 'Memory used, in MB:'
   write(file_id,222) 'cell stencils:                ',size(cellstencils)/1024.**2    
   write(file_id,222) 'meshelements:                 ',size(meshelements)/1024.**2
   write(file_id,222) 'solution variable CSV and PV: ',(size(csv)+size(csv))/1024.**2
   write(file_id,222) 'fluxes:                       ',size(flux)/1024.**2
   write(file_id,222) 'degrees of freedom:           ',size(DegreesOfFreedom)/1024.**2
   write(file_id,222) 'extrapolateddegreesOffreedom: ',size(ExtrapolatedDegreesOfFreedom)/1024.**2   
   write(file_id,222) 'ExtValues:                    ',size(ExtValues)/1024.**2      
   write(file_id,222) 'ExtVertexValues:              ',size(ExtVertexValues)/1024.**2          
   write(file_id,222) 'Total memory, mbs:            ',totalmemory/1024.**2
   write(file_id,*) '    '
   
   write(file_id,222) 'Memory used relative to the total memory:'   
   write(file_id,222) 'cell stencils :                ',size(cellstencils)/totalmemory    
   write(file_id,222) 'meshelements :                 ',size(meshelements)/totalmemory
   write(file_id,222) 'solution variable CSV and PV : ',(size(csv)+size(csv))/totalmemory
   write(file_id,222) 'fluxes :                       ',size(flux)/totalmemory
   write(file_id,222) 'degrees of freedom :           ',size(DegreesOfFreedom)/totalmemory
   write(file_id,222) 'ExtrapolatedDegreesOfFreedom : ',size(ExtrapolatedDegreesOfFreedom)/totalmemory   
   write(file_id,222) 'ExtValues :                    ',size(ExtValues)/totalmemory      
   write(file_id,222) 'ExtVertexValues :              ',size(ExtVertexValues)/totalmemory          

   write(file_id,*)'-------------'
   write(file_id,'(2x,a,1x,e11.4,a)') 'Current time :              ', t_/time,' %'
   write(file_id,'(2x,a,1x,i)') 'Total number of time steps: :              ',it
   close(file_id)


  End subroutine

 !%%%%%%%%%%%%%%%% Calculation of the sound speed  on the conserved vector CDS
  Real function ComputeSoundSpeed(cds)
    Real cds(3),p,u  ,GM
    GM=1.4
    u = cds(2)/cds(1)
    p = (GM-1)*(cds(3) - 0.5*cds(2)*u)
    ComputeSoundSpeed=sqrt(gm*p/cds(1))  
  End function



!  compute the Godunov flux with the HLLC Riemann solver
   SUBROUTINE HLLC_FLUX(CDL,CDR,FHLLC)

      INTEGER  K 
      Real(8), DIMENSION(4) :: FHLLC(4)
      Real(8)  ENEL, ENER,SL,SM,SR
      Real(8), DIMENSION(4) ::  FDL, FDR, FSL, FSR
      Real(8), DIMENSION(4) ::  CSL, CSR, CDL, CDR,PVL,PVR
      REAL(8) DL,UL,VL,PL,DR,UR,VR,PR,CL,CR,EL,ER

      Real  ULstar(4),URstar(4),FL(4), FR(4),Flux(4),AL,AR,S_star
	
      PVL = csvtopv(CDL)
      PVR = csvtopv(CDR)
	CALL ESTIME(PVL,PVR,SL,SR,SM)
        FDL = FLUEVAL(CDL)
        FDR = FLUEVAL(CDR) 
	      DL = CDL(1)
	      DR = CDR(1)
         UL = CDL(2)/CDL(1)
         UR = CDR(2)/CDR(1)
	       VL = CDL(3)/CDL(1)
	       VR = CDR(3)/CDR(1)
 	       PL = FDL(2)-(CDL(2)*(CDL(2)/CDL(1)))
         PR = FDR(2)-(CDR(2)*(CDR(2)/CDR(1)))
	       EL = CDL(4) 
         ER = CDR(4) 
         

        S_star= ((PR-PL)+((DL*(UL))*(SL-(UL)))-((DR*(UR))*(SR-(UR))))/(((DL)*(SL-(UL)))-(DR*(SR-(UR))))

         ULstar(1)=DL*((SL-(UL))/(SL-S_star))
         ULstar(2)=DL*((SL-(UL))/(SL-S_star))*S_star
         ULstar(3)=DL*((SL-(UL))/(SL-S_star))*VL
	       ULstar(4)=DL*((SL-(UL))/(SL-S_star))*((EL/DL)+((S_star-(UL))*(S_star+(PL/(DL*(SL-(UL)))))))
	
         URstar(1)=DR*((SR-(UR))/(SR-S_star))
         URstar(2)=DR*((SR-(UR))/(SR-S_star))*S_star
         URstar(3)=DR*((SR-(UR))/(SR-S_star))*VR
	       URstar(4)=DR*((SR-(UR))/(SR-S_star))*((ER/DR)+((S_star-(UR))*(S_star+(PR/(DR*(SR-(UR)))))))
	

         if (SL.GE.0) then 
            FHLLC=FDL
         else if (SR.LE.0)then 
            FHLLC=FDR
         else if (SL.LE.0.and.S_star.GE.0) then 
            FHLLC=FDL+(SL*(ULstar-CDL))
         else if (S_star.LE.0.and. SR.GE.0) then 
            FHLLC=FDR+(SR*(URstar-CDR))
         end if 

     ! Estimates for wave speeds are already calculated

     
	
   END SUBROUTINE


!  compute the Godunov flux with the HLL Riemann solver
   SUBROUTINE RUSANOV_FLUX(CDL,CDR,FHLL)

      INTEGER  K 
      Real(8), DIMENSION(4) :: FHLL(4)
      Real(8)  ENEL, ENER,SL,SR,SM
      Real(8), DIMENSION(4) ::  FDL, FDR, FSL, FSR
      Real(8), DIMENSION(4) ::  CSL, CSR, CDL, CDR,PVL,PVR
  
      PVL = csvtopv(CDL)
      PVR = csvtopv(CDR)

      CALL ESTIME(PVL,PVR,SL,SR,SM)

     ! Estimates for wave speeds are already calculated

!      SL = min(PVL(2)-dsqrt(gamma*PVL(4)/PVL(1)),PVR(2)-dsqrt(gamma*PVR(4)/PVR(1)))
!      SR = max(PVL(2)+dsqrt(gamma*PVL(4)/PVL(1)),PVR(2)+dsqrt(gamma*PVR(4)/PVR(1)))      

     ! fully supersonic flow to the left
      	 
       ! Rusanov flux
        FDL = FLUEVAL(CDL)
        FDR = FLUEVAL(CDR)        
        FHLL = 0.5*(FDL+FDR) - 0.5*max(abs(SL),abs(SR))*(CDR-CDL)
!       return
	 
   END SUBROUTINE



! *************** Construction of stenicls *********************************************
  Subroutine PrestoreReconstructionData
    Integer I,k,L,j,kand,m,iii,transl
    Integer StencilNum,StencilSizeCounter
    Real(8) xco(2),xs(2),XA(2),XB(2),XC(2),v1(2),v2(2),v3(2),v4(2),S,test
    Real(8) Jacobian(2,2), InverseJacobian(2,2),Determ,Sj
    Integer(4), allocatable :: CandStencil(:) ! local values used during the search
    Integer CandStencilSize, NumOfCand
    Type(CandidateCell) ListOfcand(30)
    
    Integer     MaxStencilSize, ActualStencilSize

    Integer(4), allocatable:: Stencil(:,:),StencilSize(:)
    Real(8), allocatable   :: LocalElementCentres(:,:,:),LocalElementVertexes(:,:,:,:),LocalElementArea(:,:)
    Real(8), allocatable   :: MA(:,:,:),  Akl(:,:,:),ExtValues(:), InvertedMA(:,:,:)
    Integer, allocatable   :: LocalElementNumSides(:,:)
    Real(8), allocatable   :: IntBasisFunction(:)
    Real(8), allocatable   :: IndicatorMatrix(:,:)

    Logical  MixedElement

    Real(8), allocatable :: Trans(:,:) 

    Integer AllCellsCounter
    Integer AllCellLists(500)

    
   print*,' Start prestore'


    iii = NElementMax/5

    ! 1. Form the stencils in physical coordinates

    Do i=1, NElementMax


     ! select the numer of cells in the stencil
     MaxStencilSize = int(max(1.D0*MeshElements(i)%NumberOfSides,StencilSizeMultiplayer*TheoNumDeg))

    ! allocate all arrays
     allocate(Stencil(0:MaxNumberOfStencils,0:MaxStencilSize),StencilSize(0:MaxNumberOfStencils))
     allocate(Akl(0:MaxNumberOfStencils,MaxStencilSize,MaxNumbDegrees), &
            LocalElementCentres(0:MaxNumberOfStencils,0:MaxStencilSize,2), &
	    LocalElementNumSides(0:MaxNumberOfStencils,0:MaxStencilSize), &
            LocalElementVertexes(0:MaxNumberOfStencils,0:MaxStencilSize,4,2), &
            LocalElementArea(0:MaxNumberOfStencils,0:MaxStencilSize))
     allocate(MA(0:MaxNumberOfStencils,MaxNumbDegrees,MaxNumbDegrees))
     allocate(InvertedMA(0:MaxNumberOfStencils,MaxNumbDegrees,MaxNumbDegrees))
     allocate(IntBasisFunction(0:MaxNumbDegrees))
     allocate(IndicatorMatrix(MaxNumbDegrees,MaxNumbDegrees))
     allocate(CandStencil(0:MaxStencilSize))


    ! ----------------------------------------------------------------------------
    ! 1a. Form the central stencil

   !  print*,' i=',i
     CellStencils(i)%iNumberOfStencils = 0

     StencilSize(:) = 0
     CandStencil(:) = 0
     CandStencil(0) = i
     StencilSizeCounter = 0 

     StencilNum = 0 ! - not to forget

     !	i) loop over all cells and create the list of neigbours
    
     5551 Continue

     NumOfCand = 0
     Do k=0,StencilSizeCounter
        Do L=1,MeshElements(CandStencil(k))%NumberOfSides	   
         ! candidate cell number
          kand =  MeshElements(CandStencil(k))%NeighborsID(L) 

          if (kand <0) then
            ! the side of the element is attached to a solid wall
            ! skip  the element
            goto 101
          endif
          
 	      ! 1. check  that it is not already in the stencil or in the candidate's list
          do j=0,StencilSizeCounter
	          If (kand .eq. CandStencil(j)) goto 101
          enddo
          do j=1, NumOfCand
	          If (kand .eq. ListOfCand(j)%id) goto 101
          enddo
           ! 2. all checks are passed  - add the cell into the stencil
!          print*,' add cell'
          NumOfCand = NumOfCand + 1
          ListOfCand(NumofCand)%id = kand
          xs = MeshElements(i)%center
          xco = MeshElements(kand)%center
          ! account for periodic boundaries 
	      if (abs(xco(1) - xs(1)) .ge. XPer/2)     xco(1) = xco(1) + XPer*sign(1.D0,xs(1) - xco(1))
  	      if (abs(xco(2) - xs(2)) .ge. YPer/2)     xco(2) = xco(2) + YPer*sign(1.D0,xs(2) - xco(2))
          ListOfCand(NumOfCand)%distance = distance(xs,xco)
          101 continue
        Enddo ! L=1,NumElementangleSides(i)
     Enddo ! k =1 

     ! ii)  check possible options
     
   !    attention here!!!
     Select case(StencilSizeCounter + NumOfCand .le.  MaxStencilSize)
      case(.true.)
         ! add all cells from the list to the stencil 
         do k=1,NumOfCand
            CandStencil(StencilSizeCounter+k) =  ListOfCand(k)%id
         enddo
         StencilSizeCounter = StencilSizeCounter + NumOfCand
      case(.false.)
         ! add only closest cells
         ! reorder the array 
         CALL BubbleSort(NumOfCand,ListOfCand(1:NumOfCand))
         ! add the minimum required number of cells to the stencil 
        ! attention here!!!
         do k=1,(MaxStencilSize-StencilSizeCounter)
            CandStencil(StencilSizeCounter+k) =  ListOfCand(k)%id
         enddo
         StencilSizeCounter = MaxStencilSize ! MaxStencilSize34(MeshElements(i)%NumberOfSides)
     End select
    
     If (StencilSizeCounter < MaxStencilSize ) then ! MaxStencilSize34(MeshElements(i)%NumberOfSides)) then
 !       print*, StencilSizeCounter,' before goto 5551'
        Goto 5551 ! carry on addition of cells
     endif

     ! finally store the found central stencil
     Stencil(StencilNum,0:StencilSizeCounter) = CandStencil(0:StencilSizeCounter)
     StencilSize(StencilNum) = StencilSizeCounter
 !    print*,' final size ',StencilSizeCounter,' of the stencil', StencilNum
     Do k=0,StencilSizeCounter
        LocalElementCentres(StencilNum,k,:)    = MeshElements(CandStencil(k))%Center 
        LocalElementNumSides(StencilNum,k)     = MeshElements(CandStencil(k))%NumberofSides 
        do m=1,MeshElements(CandStencil(k))%NumberOfSides 
            LocalElementVertexes(StencilNum,k,m,:) =  MeshElements(CandStencil(k))%VertexCoordinates(:,m)    
        enddo ! m
        LocalElementArea(StencilNum,k)  = MeshElements(CandStencil(k))%Area 
     Enddo

    ! ----------------------------------------------------------------------------
    ! 1b. Form the sectorial stencils
    !     we use the sectors, formed by two successive verteces and the cell center
    
     AllCellsCounter = 0
     AllCellLists = 0


   ! if the scheme is not WENO, skip this part
       If (SchemeType <3) Goto 4444

!     print*,'Form the sectorial stencils'
     
     CellStencils(i)%iNumberOfStencils = 0
     
    ! Loop over all sectorial stencils
     Do StencilNum =1,MeshElements(i)%NumberOfSides

      CandStencil(:) = 0
      CandStencil(0) = i
      StencilSizeCounter = 0 

     ! sectorial stencils - precompute transformation matrixes for searching
      xs = MeshElements(i)%Center ! center of the main cell i
      Select case(MeshElements(i)%NumberOfSides)
	       case(3)
	       ! triangular cell
            v1  = xs
            v2  = MeshElements(i)%VertexCoordinates(:,VertexToSide3(StencilNum,2))
            v3  = MeshElements(i)%VertexCoordinates(:,VertexToSide3(StencilNum,3))
	       case(4)
	        ! quadrilateral cell
           v1  = xs
           v2  = MeshElements(i)%VertexCoordinates(:,VertexToSide4(StencilNum,2))
           v3  = MeshElements(i)%VertexCoordinates(:,VertexToSide4(StencilNum,3))
	  End select ! Stencil Number
	
      CALL ComputeJacobians(v1,v2,v3,Jacobian,InverseJacobian,Determ)

     !	i) loop over all cells and create the list of neigbours
      5552 Continue
      NumOfCand = 0
      Do k=0,StencilSizeCounter
        Do L=1,MeshElements(CandStencil(k))%NumberOfSides	   
         ! candidate cell number
          kand =  MeshElements(CandStencil(k))%NeighborsID(L) 
          
          if (kand <0) then
            ! the side of the element is attached to a solid wall
            ! skip  the element
            goto 201
          endif

          
 	      ! 1. check  that it is not already in the stencil or in the candidate's list
          do j=0,StencilSizeCounter
	          If (kand .eq. CandStencil(j)) goto 201
          enddo
          do j=1, NumOfCand
	          If (kand .eq. ListOfCand(j)%id) goto 201
          enddo

	      ! 2. check that cell centre lies within the sector; for sectorial stencils only
          XC = MeshElements(kand)%Center
             ! take into account periodic boundaries
          if (abs(xc(1) - xs(1)) .ge. XPer/2)    xc(1) = xc(1) + XPer*sign(1.D0,xs(1) - xc(1))
          if (abs(xc(2) - xs(2)) .ge. YPer/2)    xc(2) = xc(2) + YPer*sign(1.D0,xs(2) - xc(2))
          XC =  matmul(InverseJacobian(:,:),XC -v1)
          If ( XC(1) .le. -0.05 .or. XC(2) .le. -0.05)   Goto 201

 	      ! 3. check  that it is not already used in some other stencils
          do j=1,AllCellsCounter
	          If (kand .eq. AllCellLists(j)) then
	             goto 201
              endif	             
          enddo

          ! 4. all checks are passed  - add the cell into the stencil
          NumOfCand = NumOfCand + 1
          ListOfCand(NumofCand)%id = kand
          xs = MeshElements(i)%center
          xco = MeshElements(kand)%center
          ! account for periodic boundaries 
	      if (abs(xco(1) - xs(1)) .ge. XPer/2)     xco(1) = xco(1) + XPer*sign(1.D0,xs(1) - xco(1))
  	      if (abs(xco(2) - xs(2)) .ge. YPer/2)     xco(2) = xco(2) + YPer*sign(1.D0,xs(2) - xco(2))
          ListOfCand(NumOfCand)%distance = distance(xs,xco)
          
          ! all to the global list of cells in secttorail stencils
          AllCellsCounter = AllCellsCounter+1
          AllCellLists(AllCellsCounter) =  kand

          201 continue
        Enddo ! L=1,NumElementangleSides(i)
      Enddo ! k =1, StMax

     ! ii)  check possible options
     
      ! first check if we can add any candidates
      If (NumOfCand .eq. 0) then
        ! it is not possible to form a sectorial stencil here
        ! go to the next side
        Goto 5553
      Endif  

      Select case(StencilSizeCounter + NumOfCand .le. MaxStencilSize) !  MaxStencilSize34(MeshElements(i)%NumberOfSides) )
       case(.true.)
         ! add all cells from the list to the stencil 
         do k=1,NumOfCand
            CandStencil(StencilSizeCounter+k) =  ListOfCand(k)%id
         enddo
         StencilSizeCounter = StencilSizeCounter + NumOfCand
       case(.false.)
         ! add only closest cells
         ! reorder the array 
         CALL BubbleSort(NumOfCand,ListOfCand(1:NumOfCand))
         ! add the minimum required number of cells to the stencil 
         do k=1,(MaxStencilSize-StencilSizeCounter)
            CandStencil(StencilSizeCounter+k) =  ListOfCand(k)%id
         enddo
         StencilSizeCounter = MaxStencilSize
      End select
    
      If (StencilSizeCounter < MaxStencilSize)  Goto 5552 ! carry on addition of cells

      ! the stencil is formed & stored

      ! increase the number of stencils corresponding to cell i
      CellStencils(i)%iNumberOfStencils = CellStencils(i)%iNumberOfStencils + 1

      Stencil(CellStencils(i)%iNumberOfStencils ,0:StencilSizeCounter) = CandStencil(0:StencilSizeCounter)
      StencilSize(CellStencils(i)%iNumberOfStencils) = StencilSizeCounter           
!      print*,' final size ',StencilSizeCounter,' of the stencil', StencilNum
      Do k=0,StencilSizeCounter
        LocalElementCentres(CellStencils(i)%iNumberOfStencils,k,:)    = MeshElements(CandStencil(k))%Center 
        LocalElementNumSides(CellStencils(i)%iNumberOfStencils,k)     = MeshElements(CandStencil(k))%NumberofSides 
        do m=1,MeshElements(CandStencil(k))%NumberOfSides 
                LocalElementVertexes(CellStencils(i)%iNumberOfStencils,k,m,:) =   &
                      MeshElements(CandStencil(k))%VertexCoordinates(:,m)    
        enddo ! m
        LocalElementArea(CellStencils(i)%iNumberOfStencils,k)  = MeshElements(CandStencil(k))%Area 
     Enddo

     5553 Continue   ! we jump to this lable if no stencil can be formed
 
    Enddo ! End of the Stencil Num loop


   4444 Continue  ! this is used to avoid using sectorial stencils for linear and Tvd of Barth schemes

     ! ************** account for periodic boundaryies
     ! ************** If a certain cell is on the other side of the boundary shift the cell coordinates


 !    print*,'account for periodic boundaryies'
     Do StencilNum = 0,CellStencils(i)%iNumberOfStencils ! loop over stencils for a given cell 'i'
  !     print*,'StencilNum=',StencilNum
  !     print*,'StencilSize(StencilNum,i) =',StencilSize(StencilNum,i)
       xs = LocalElementCentres(StencilNum,0,:)
       Do k=1,StencilSize(StencilNum)
!           print*,' k=',k
        ! check if the cell centre positions in the stencil are close to the boundary
	     xco = LocalElementCentres(StencilNum,k,:)
	     if (abs(xco(1) - xs(1)) .ge. XPer/2) then
             transl = sign(1.D0,xs(1) - xco(1))
  	         xco(1) = xco(1) + XPer*transl
      	     do L=1,LocalElementNumSides(StencilNum,k)
	              LocalElementVertexes(StencilNum,k,L,1) = &
                  LocalElementVertexes(StencilNum,k,L,1) + XPer*transl
             enddo
          endif

	     if (abs(xco(2) - xs(2)) .ge. YPer/2)  then
             transl = sign(1.D0,xs(2) - xco(2))
             xco(2) = xco(2) + YPer*transl
             do L=1,LocalElementNumSides(StencilNum,k)
                LocalElementVertexes(StencilNum,k,L,2) = &
                       LocalElementVertexes(StencilNum,k,L,2) + YPer*transl
             enddo
         endif
         LocalElementCentres(StencilNum,k,:) = xco
       Enddo

     Enddo ! stenciulnum
    

   ! ************  Now transform all stencils to the reference state ************

  !  print*,'transform all stencils to the reference state'
      ! define the transformation for cell i
     Select case(MeshElements(i)%NumberOfSides)
         case(3)
          v1  = LocalElementVertexes(0,0,1,:)
          v2  = LocalElementVertexes(0,0,2,:)
          v3  = LocalElementVertexes(0,0,3,:)
	 case(4)
          v1  = LocalElementVertexes(0,0,1,:)
          v2  = LocalElementVertexes(0,0,2,:)
          v3  = LocalElementVertexes(0,0,4,:)
      End select
  
     CALL ComputeJacobians(v1,v2,v3,Jacobian,InverseJacobian,Determ)

     Do StencilNum = 0,CellStencils(i)%iNumberOfStencils ! loop over stencils for a given cell 'i'
      ! Transform all coordinates to the reference state
       Do k=0,StencilSize(StencilNum)
          LocalElementCentres(StencilNum,k,:) = matmul(InverseJacobian(:,:),LocalElementCentres(StencilNum,k,:)-v1)
          Do L=1,MeshElements(Stencil(StencilNum,k))%NumberOfSides
           LocalElementVertexes(StencilNum,k,L,:) = &
                 Matmul(InverseJacobian(:,:),LocalElementVertexes(StencilNum,k,L,:) - v1)
          Enddo
          LocalElementArea(StencilNum,k) =LocalElementArea(StencilNum,k)/Determ
       Enddo
     Enddo ! End of StencilNum loop



   ! stencil, all coordinates & areas are in reference coordinate system

  !   print*,' store stencil, all coordinates & areas are in reference coordinate system'

!    CellStencils(i)%iNumberOfStencils =  MeshElements(i)%NumberOfStencils
    allocate(CellStencils(i)%Stencils(0:CellStencils(i)%iNumberOfStencils))

    CellStencils(i)%iNumberOfSides      = MeshElements(i)%NumberOfSides
    CellStencils(i)%iNumberOfStencils   = CellStencils(i)%iNumberOfStencils

    allocate(CellStencils(i)%rCellVertexes(CellStencils(i)%iNumberOfSides,2))
    allocate(CellStencils(i)%rIntBasisFunctions(0:MaxNumbDegrees))
    allocate(CellStencils(i)%rIndicatorMatrix(MaxNumbDegrees,MaxNumbDegrees))

    CellStencils(i)%rCellVertexes(:,:) = LocalElementVertexes(0,0,:,:)
    CellStencils(i)%rIntBasisFunctions = 0.D0
    CellStencils(i)%rIndicatorMatrix   = 0.D0

    CellStencils(i)%rJacobian  = Jacobian
    CellStencils(i)%rInverseJacobian  = InverseJacobian

    do StencilNum=0,CellStencils(i)%iNumberOfStencils

    ! here check the types of cells in the stencil
    ! if all cells are of the same type than decrease the size of the stencil from 2 to 1.5
 
    MixedElement = .false.
    Do k=1,MaxStencilSize
     L = Stencil(StencilNum,k)
     if (MeshElements(L)%NumberOfSides .ne.  MeshElements(Stencil(StencilNum,0))%NumberOfSides)  MixedElement = .true.
    Enddo

    Select case(MixedElement)
      case(.true.)
!        print*,' mixed element detected for i=',i
       ActualStencilSize =  StencilSize(StencilNum)
      case(.false.)
!       print*,' StencilSize(StencilNum) = ',StencilSize(StencilNum)
!      print*,' reduce the stencil from ',MaxStencilSize, &
!          ' elements to ',  int(max(1.D0*MeshElements(i)%NumberOfSides,1.5*TheoNumDeg))
!       read*
       ActualStencilSize = int(max(1.D0*MeshElements(i)%NumberOfSides,1.5*TheoNumDeg))      
    End select

    ! maximum stencil size for array allocation
!    StencilSizeMultiplayer = 2.
     ! select the numer of cells in the stencil
!     MaxStencilSize = int(max(1.D0*MeshElements(i)%NumberOfSides,StencilSizeMultiplayer*TheoNumDeg))

     CellStencils(i)%Stencils(StencilNum)%iSize  = ActualStencilSize 

     allocate(CellStencils(i)%Stencils(StencilNum)%iCellIDs(0:CellStencils(i)%Stencils(StencilNum)%iSize))
     allocate(CellStencils(i)%Stencils(StencilNum)%rElementCentres(0:CellStencils(i)%Stencils(StencilNum)%iSize,2))
     allocate(CellStencils(i)%Stencils(StencilNum)%rElementVertexes(0:CellStencils(i)%Stencils(StencilNum)%iSize,4,2))
     allocate(CellStencils(i)%Stencils(StencilNum)%rElementAreas(0:CellStencils(i)%Stencils(StencilNum)%iSize))
     allocate(CellStencils(i)%Stencils(StencilNum)%rAkl(CellStencils(i)%Stencils(StencilNum)%iSize,MaxNumbDegrees))
     allocate(CellStencils(i)%Stencils(StencilNum)%iElementNumSides(0:CellStencils(i)%Stencils(StencilNum)%iSize))


     allocate(CellStencils(i)%Stencils(StencilNum)%rMA(MaxNumbDegrees,MaxNumbDegrees))
     CellStencils(i)%Stencils(StencilNum)%rMA = 0.

     allocate(CellStencils(i)%Stencils(StencilNum)%rQmatrix(MaxNumbDegrees,MaxNumbDegrees))
     CellStencils(i)%Stencils(StencilNum)%rQmatrix = 0.

     allocate(CellStencils(i)%Stencils(StencilNum)%rRmatrix(MaxNumbDegrees,MaxNumbDegrees))
     CellStencils(i)%Stencils(StencilNum)%rRmatrix = 0.

     allocate(CellStencils(i)%Stencils(StencilNum)%rInvMatrix(MaxNumbDegrees,MaxNumbDegrees))
     CellStencils(i)%Stencils(StencilNum)%rInvMatrix = 0.

     allocate(CellStencils(i)%Stencils(StencilNum)%rFinalMatrix(MaxNumbDegrees,CellStencils(i)%Stencils(StencilNum)%iSize))
     CellStencils(i)%Stencils(StencilNum)%rFinalMatrix = 0.


     CellStencils(i)%Stencils(StencilNum)%iCellIDs(0:CellStencils(i)%Stencils(StencilNum)%iSize) = &
              Stencil(StencilNum,0:CellStencils(i)%Stencils(StencilNum)%iSize)

     CellStencils(i)%Stencils(StencilNum)%rElementCentres(0:CellStencils(i)%Stencils(StencilNum)%iSize,:) = &
            LocalElementCentres(StencilNum,0:CellStencils(i)%Stencils(StencilNum)%iSize,:)

     CellStencils(i)%Stencils(StencilNum)%rElementVertexes(0:CellStencils(i)%Stencils(StencilNum)%iSize,:,:) = &
          LocalElementVertexes(StencilNum,0:CellStencils(i)%Stencils(StencilNum)%iSize,:,:)

     CellStencils(i)%Stencils(StencilNum)%rElementAreas(0:CellStencils(i)%Stencils(StencilNum)%iSize) = &
            LocalElementArea(StencilNum,0:CellStencils(i)%Stencils(StencilNum)%iSize)

     CellStencils(i)%Stencils(StencilNum)%iElementNumSides(0:CellStencils(i)%Stencils(StencilNum)%iSize) = &
          LocalElementNumSides(StencilNum,0:CellStencils(i)%Stencils(StencilNum)%iSize)
   Enddo ! StencilNum

   ! since local coordinates are defined, compute integrals
    CellStencils(i)%rArea = LocalElementArea(0,0)
    CellStencils(i)%rIntBasisFunctions = CalcIntBasisFunction(i)

  ! compute indicator
    If (schemetype   .eq. 3) then
!         print*,' schemetype=',schemetype 
         CALL ComputeIndicatorMatrix(i,CellStencils(i)%rIndicatorMatrix(:,:))
    else 
         CellStencils(i)%rIndicatorMatrix(:,:) = 0.D0
    endif


   ! deallocate arrays
  deallocate(Stencil,StencilSize)
  deallocate(Akl, LocalElementCentres, LocalElementNumSides,  LocalElementVertexes, &
            LocalElementArea, MA, InvertedMA,IntBasisFunction,IndicatorMatrix)
  deallocate(CandStencil)

  Enddo ! i=1,NelementMax


   print*,' precompute integrals for LSQ reconstruction'
	

    !******************** precompute integrals for LSQ reconstruction ******************

    ! this part is called for high order schemes only

    If (SpatialOrder > 1) then 


    iii = NElementMax/5
    write(*,'(2x,a)')' progress: '

    Do i=1,NElementMax

     if ( mod(i,iii) .eq. 0 .and. i>1) write(*,'(2x,f6.2)')1.*i/NElementMax

     Do StencilNum = 0,CellStencils(i)%iNumberOfStencils

       CellStencils(i)%Stencils(StencilNum)%rAkl = 0.
       
       Do k=1,MaxNumbDegrees

        ! the integral of the basis function of the element is substructed automatically
		! calculate Akl, which is spatial average of the basis function "k"
		! over the local element "j"
         Do j=1,CellStencils(i)%Stencils(StencilNum)%iSize 

          ! coordinates of vertexes of the cell "j"
         select case(CellStencils(i)%Stencils(StencilNum)%iElementNumSides(j))
          case(3)
           v1 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,1,:)
           v2 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,2,:)
           v3 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,3,:)
           case(4)
           v1 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,1,:)
           v2 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,2,:)
           v3 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,3,:)
           v4 = CellStencils(i)%Stencils(StencilNum)%rElementVertexes(j,4,:)
          end select

           xa = CellStencils(i)%Stencils(StencilNum)%rElementCentres(0,:)
           xb = CellStencils(i)%Stencils(StencilNum)%rElementCentres(j,:)

		   ! spatial integral
          CellStencils(i)%Stencils(StencilNum)%rAkl(j,k)  = &
     ComputeBasisFunctionIntegralsOverElement(v1,v2,v3,v4,CellStencils(i)%Stencils(StencilNum)%iElementNumSides(j),k,i)

	       ! divide over cell size to obtain the cell average
		   sj = CellStencils(i)%Stencils(StencilNum)%rElementAreas(j)
     	   CellStencils(i)%Stencils(StencilNum)%rAkl(j,k)  = &
		      CellStencils(i)%Stencils(StencilNum)%rAkl(j,k)/sj
	  Enddo
      Enddo ! loop over degrees of freedom

       !   form the  coefficient matrix for least squares Ma
	   ! element k,l of the matrix is calculated as the sum
	   ! Ma(k,l) = sum_(j=1,size of the stencil) Akl(j,k) Akl(j,l)
        
       CellStencils(i)%Stencils(StencilNum)%rMA = 0.
       do L=1,MaxNumbDegrees
           do k=1,MaxNumbDegrees
              ! loop over all stencils
              Do j=1,CellStencils(i)%Stencils(StencilNum)%iSize 
                  xa = CellStencils(i)%Stencils(StencilNum)%rElementCentres(0,:)
                  xb = CellStencils(i)%Stencils(StencilNum)%rElementCentres(j,:)

                  ! here we need to divide over SJ otherwise it is used twice
                  CellStencils(i)%Stencils(StencilNum)%rMA(k,L) = CellStencils(i)%Stencils(StencilNum)%rMA(k,L) + &
           CellStencils(i)%Stencils(StencilNum)%rAkl(j,k)*CellStencils(i)%Stencils(StencilNum)%rAkl(j,L)
               enddo          
            enddo
       Enddo

      ! perform the QR decomposition of the matrix Ma
	  ! and also invert the QR matrix
	  ! rInvMatrix = (QR)^(-1)
      CALL QRDecomposition(CellStencils(i)%Stencils(StencilNum)%rMA, & 
                           CellStencils(i)%Stencils(StencilNum)%rQmatrix, & 
                           CellStencils(i)%Stencils(StencilNum)%rRmatrix, &
                           CellStencils(i)%Stencils(StencilNum)%rInvMatrix, &
                           MaxNumbDegrees)

       !    finalsolution        
		allocate(Trans(MaxNumbDegrees,CellStencils(i)%Stencils(StencilNum)%iSize))

		! transpose the matrix Akl of integrals of basis functions over cells in the stencil 
        do L=1,MaxNumbDegrees
           do j=1,CellStencils(i)%Stencils(StencilNum)%iSize
		     Trans(L,j) =CellStencils(i)%Stencils(StencilNum)%rAkl(j,L)
           enddo
        enddo		

		! calculate the final matrix, used later for evaluation of the degree of freedom "dw"
		! the least square system is Ma*dw = Ark^T (linear function of cell averages)
		! or
		!  Q*R*dw = Trans*()
		! or
		! dt = Final *(), where Final = (QR)^{-1)*Trans

        CellStencils(i)%Stencils(StencilNum)%rFinalMatrix = &
		  matmul(CellStencils(i)%Stencils(StencilNum)%rInvmatrix(:,:),Trans)

		deallocate(Trans)
  

     Enddo ! stencilnum

   Enddo ! i=1,NElementMax

   Endif ! spatail order > 1

   ! write out stencils in reference space
!   Call StencilOutput(1)
!   Call StencilOutput(2)

   print*,' Done prestore'
   print*

  End Subroutine

  Subroutine WriteRestart
    integer file_id,i
    
   print*,' writing into restart file'
   file_id = 555
   open(file_id,file='restart.bin')
   do i=1,NELementMax
     write(file_id,*)csv(1,:,i)
   enddo
   write(file_id,*)t_
   write(file_id,*)it
   close(file_id)

   print*,' writing into dublicate file file'
   file_id = 555
   open(file_id,file='restart2.bin')
    do i=1,NELementMax
     write(file_id,*)csv(1,:,i)
    enddo
    write(file_id,*)t_
    write(file_id,*)it
   close(file_id)

   print*,' writing complete'
      
  End Subroutine



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Function InitialData(x,y)
    Real(8) x,y,r, InitialData(4),eps,n,u,v,T,p, ene,PV0(4)
    Real(8) ::Dv=1, V0 = 0.1, Minf=0.1
    Real(8), parameter :: pi= 3.141592653589793, kh = 2*pi   
    Real(8)     DenPostShock, UPostShock,PPostShock,ShockSpeed
    Real(8) ar,mr,s1,S2
    Real(8) p1,p2
    Real(8), parameter :: DenR = 1.4, VelR = 0., PresR=1., MS=10. 
     
!     InitialV(1) = 1.
!     InitialV(2) = 2.
!     InitialV(3) = -3.	 
!     InitialV(4) = 1.
!     return

 !   print*,' IC:',inicondtype    read*


   Select case(IniCondType)
    case(1) ! vortex     

     if (x .le. xmin)  x = x + xper
     if (y .le. ymin)  y = y + yper
     if (x .ge. xmax)  x = x - xper
     if (y .ge. ymax)  y = y - yper

     r = dsqrt(x**2+y**2)

     eps = 5.
     u = 1 - (eps/2/pi)*exp(0.5*(1-r**2))*y
     v = 1 + (eps/2/pi)*exp(0.5*(1-r**2))*x
     T = 1 - (gamma-1)*(eps**2/8/pi**2/gamma)*exp(1-r**2)
     n = T**(1./(gamma-1))
     p = n*T

    case(2) ! explosion

     u = 0.D0   
     v = 0.D0
     r = dsqrt(x**2+y**2)
  !   r = abs(x)
     if (r .le. 0.4) then
        n = 1.D0
        p = 1.D0
     else
        n = 0.125D0;    p = 0.1D0    
     endif
     
     
     case(3) ! implosion

     u = 0.D0   
     v = 0.D0
     r = dsqrt(x**2+y**2)
  !   r = abs(x)
     if (r .gt. 0.4) then
        n = 1.D0
        p = 1.D0
     else
        n = 0.125D0;    p = 0.1D0    
     endif
     
     
     
     case(5) ! forward facing step

     u = 3.D0   
     v = 0.D0
     n = 1.0d0
  !   r = abs(x)
     p=1.0d0/1.4d0
     
     case(6) ! double mach reflection
     
     if (x.lt.((1.0d0/6.0d0)+(y/(dsqrt(3.0d0)))))then
        n=8.0d0
        u=8.25*cos(pi/6.0d0)
        v=-8.25*sin(pi/6.0d0)
        p=116.5
        
    else
        n=1.4d0
        u=0.0D0
        v=0.0D0
        p=1.0d0
        
    end if
    
    
     case(0)

     n = 1.;     u = 0.;     v = 0.;     p = 1.


!      case(3)
! 
! !     if (x<1./6) then
!      if (x<-0.1     ) then     
! !       AR = sqrt(gamma*PresR/DenR)
! !       MR = VelR/Ar
! !       PPostShock   = PresR*( 2*gamma*(Mr-Ms)**2 - (gamma-1))/(gamma+1)
! !       S1 = Ppostshock/PresR
! !       S2 = (gamma-1)/(gamma+1)
! !       DenPostShock =  DenR*( S1 + S2 )/( S1*S2 + 1)
! !       DenPostShock = DenR*(gamma+1)*(Mr-Ms)**2/( (gamma-1)*(Mr-MS)**2 + 2)
! !       S1 = (gamma+1)/(2*gamma) 
! !       S2 = (gamma-1)/(2*gamma) 
! !       ShockSpeed   = Velr + Ar*sqrt(  S1*(Ppostshock/PresR) + S2)
! !       Upostshock = (1 - DenR/Denpostshock)*Shockspeed + VelR*DenR/DenPostShock
! !       open(444,file='postshock.txt')
! !       write(444,*) denr,0.,presr
! !       write(444,*) 'density:', denpostshock
! !       write(444,*) 'velocity:',upostshock
! !       write(444,*) 'pressure:',ppostshock
! !       write(444,*) 'temperature:',  ppostshock/denpostshock
! !       write(444,*) 'speed:',shockspeed
! !       close(444)
! 
!        DenPostShock = 8.
!        UpostShock = 8.25
!        Ppostshock = 116.5
!        ShockSpeed = 10
!        n = DenPostShock     
!        u = Upostshock  
!        v = 0.
!        p = PPostShock
!      else
!        n= DenR; u=0.; v=0.; p=PresR
!      endif
! 
! 
! !       n= 1.; u=1.; v=0.; p=1.

     case(4)

     n = 1.4;     u = 3.;     v = 0.;     p = 1.


    End select


     PV0(1) = n
     PV0(2) = u
     PV0(3) = v	 
     PV0(4) = p    
     InitialData = PvtoCSV(PV0)
  End function
 
  End
	
  
	
	
