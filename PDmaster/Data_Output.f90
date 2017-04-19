!***************************************************************
!                  Date_output                                 *
!          功能： 控制变量输出                                 *
!***************************************************************    
subroutine data_output()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    
    !if (mod(istep,mout(iload)%frequency)==0)then
    !    ii=ii+1
    !    call data_output_uniformity()
    !    maxdamage=maxval(damage)
    !    write(6,'(1x,I3,2x,I2,2x,I3,2x,f6.3,2x,f6.3)') ii,iload,iincs,loadg(iload,iincs),maxdamage
    !elseif(check_c==1)then
    !    ii=ii+1
    !    call data_output_uniformity()
    !    maxdamage=maxval(damage)
    !    write(6,'(1x,I3,2x,I2,2x,I3,2x,f6.3,2x,f6.3)') ii,iload,iincs,loadg(iload,iincs),maxdamage
    !end if
    
    maxdamage=maxval(damage)
    if (maxdamage>0.3) solve(iload)%loadtype=0
    if (maxdamage>0.3.and.mod(istep,mout(iload)%frequency)==0)then
        ii=ii+1
        call data_output_uniformity()
        write(6,'(1x,I3,2x,I2,2x,I3,2x,f6.3,2x,I4,2x,f6.3)') ii,iload,iincs,loadg(iload,iincs),istep,maxdamage
    elseif(check_c==1.OR.istep==solve(iload)%maxstep)then
        ii=ii+1
        call data_output_uniformity()
        write(6,'(1x,I3,2x,I2,2x,I3,2x,f6.3,2x,I4,2x,f6.3)') ii,iload,iincs,loadg(iload,iincs),istep,maxdamage
    end if
    
end subroutine data_output   
    
    
    
subroutine data_output_uniformity()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    call data_output_coordinates()
    call data_output_damage()
    call data_output_Displacement()
    call data_output_Velocity()
    call data_output_Force()
    call data_output_Bforce()
    call data_output_Pforce()
    !call data_output_Cforce()
    !call data_output_Gforce()
    call data_output_numfam()
    call data_output_wvolume()
    call data_output_S_fac()
    !call data_output_Strain_Energy()
    
    call data_output_case()
end subroutine data_output_uniformity
!--------------!
!   Ensight    !
!--------------!
subroutine data_output_case()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    integer i
    
    Open( 10 , File = '0point.case'  )
    write(10,'(a)')"FORMAT"
    write(10,'(a)')"type:  ensight gold"
    
    write(10,'(a)')"GEOMETRY"
    write(10,'(a)')"model:                     point.geo*       change_coords_only"
   
    write(10,'(a)')"VARIABLE"
    write(10,'(a)')"scalar per node:           Damage point.Damage*"
    write(10,'(a)')"scalar per node:           Displacement-X point.Displacement-X*"
    write(10,'(a)')"scalar per node:           Displacement-Y point.Displacement-Y*"    
	write(10,'(a)')"scalar per node:           Velocity-X point.Velocity-X*"
    write(10,'(a)')"scalar per node:           Velocity-Y point.Velocity-Y*"
    write(10,'(a)')"scalar per node:           Force-X point.Force-X*"
    write(10,'(a)')"scalar per node:           Force-Y point.Force-Y*"
    write(10,'(a)')"scalar per node:           Bforce-X point.Bforce-X*"
    write(10,'(a)')"scalar per node:           Bforce-Y point.Bforce-Y*"
    write(10,'(a)')"scalar per node:           Pforce-X point.Pforce-X*"
    write(10,'(a)')"scalar per node:           Pforce-Y point.Pforce-Y*"
    !write(10,'(a)')"scalar per node:           Cforce-X point.Cforce-X*"
    !write(10,'(a)')"scalar per node:           Cforce-Y point.Cforce-Y*"
    !write(10,'(a)')"scalar per node:           Gforce-X point.Gforce-X*"
    !write(10,'(a)')"scalar per node:           Gforce-Y point.Gforce-Y*"
    write(10,'(a)')"scalar per node:           Numfam point.Numfam*"
    write(10,'(a)')"scalar per node:           Wvolume point.Wvolume*"
    write(10,'(a)')"scalar per node:           S_fac point.S_fac*"
    !write(10,'(a)')"scalar per node:           Strain_Energy point.Strain_Energy*"
    
    write(10,'(a)')"TIME"
    write(10,'(a)')"time set:                  1"
    write(10,'(a)')"filename start number:     0"
    write(10,'(a)')"filename increment:        1"
    write(10,'(a,i4)')"number of steps:        ",ii+1
    write(10,'(a)')"time values:"
    do i=1,ii+1
        write(10,'(i6)') i
    end do
    close(10)
end subroutine data_output_case    
    
!--------------!
!     坐标     !
!--------------!
subroutine data_output_coordinates()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    real(8)::a
    integer i
    write( cStr , '(i5)' ) ii
    Open( 10 , File = 'point.geo' // AdjustL(Trim( cStr ) ) )
    write(10,'(a)')"enSight Model Geometry File"
    write(10,'(a)')"enSight 7.4.1"
    write(10,'(a)')"node id off"
    write(10,'(a)')"element id off"
    write(10,'(a)')"extents"
    write(10,'(a)')"-1.00000e+03 1.00000e+03"
    write(10,'(a)')"-1.00000e+03 1.00000e+03"
    write(10,'(a)')"-1.00000e+03 1.00000e+03"    
	write(10,'(a)')"part"
    write(10,'(a)')"         1"
    write(10,'(a)')"point"
    write(10,'(a)')"coordinates"
    write(10,'(I10)') totout(iload)
    do idime=1,3
	    do i=1,totout(iload)
           ni=outnode(iload,i)
           a=coord(ni,idime)+disp(ni,idime)     
           if(a.ge.0) then
              cFmt="(1x,e11.5)"
            else
              cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
    end do
    close(10)
end subroutine data_output_coordinates
!------------!
!近场物质点数!
!------------!    
subroutine data_output_Numfam()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    real(8)::a
    integer i
   
    write( cStr , '(i5)' ) ii
    open( 10 , File = 'point.Numfam' // AdjustL(Trim( cStr ) ) )
    write(10,'(a)')"Numfam"
    write(10,'(a)')"part"
    write(10,'(a)')"         1"
    write(10,'(a)')"coordinates"  
	do i=1,totout(iload)
        ni=outnode(iload,i)
        a=numfam(ni)    
        if(a.ge.0) then
            cFmt="(1x,e11.5)"
        else
            cFmt="(e12.5)"
        end if
        write( 10 , cFmt ) a
    end do 
    close(10)
    end subroutine data_output_Numfam
 !------------!
!近场物质点数!
!------------!    
    subroutine data_output_Wvolume()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    real(8)::a
    integer i
   
    write( cStr , '(i5)' ) ii
    open( 10 , File = 'point.Wvolume' // AdjustL(Trim( cStr ) ) )
    write(10,'(a)')"Wvolume"
    write(10,'(a)')"part"
    write(10,'(a)')"         1"
    write(10,'(a)')"coordinates"  
	do i=1,totout(iload)
        ni=outnode(iload,i)
        a=wvolume(ni)    
        if(a.ge.0) then
            cFmt="(1x,e11.5)"
        else
            cFmt="(e12.5)"
        end if
        write( 10 , cFmt ) a
    end do 
    close(10)
end subroutine data_output_Wvolume
    
!------------!
! 应变能密度 !
!------------!    
subroutine data_output_Strain_Energy()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    real(8)::a
    integer i
   
    write( cStr , '(i5)' ) ii
    open( 10 , File = 'point.Strain_Energy' // AdjustL(Trim( cStr ) ) )
    write(10,'(a)')"Strain_Energy"
    write(10,'(a)')"part"
    write(10,'(a)')"         1"
    write(10,'(a)')"coordinates"  
	do i=1,totout(iload)
        ni=outnode(iload,i)
        a=Strain_Energy(ni)    
        if(a.ge.0) then
            cFmt="(1x,e11.5)"
        else
            cFmt="(e12.5)"
        end if
        write( 10 , cFmt ) a
    end do 
    close(10)
    end subroutine data_output_Strain_Energy
    
!------------!
!表面修正系数!
!------------!    
subroutine data_output_S_fac()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    real(8)::a
    integer i
   
    write( cStr , '(i5)' ) ii
    open( 10 , File = 'point.S_fac' // AdjustL(Trim( cStr ) ) )
    write(10,'(a)')"S_fac"
    write(10,'(a)')"part"
    write(10,'(a)')"         1"
    write(10,'(a)')"coordinates"  
	do i=1,totout(iload)
        ni=outnode(iload,i)
        a=S_fac(ni)    
        if(a.ge.0) then
            cFmt="(1x,e11.5)"
        else
            cFmt="(e12.5)"
        end if
        write( 10 , cFmt ) a
    end do 
    close(10)
    end subroutine data_output_S_fac   
    
    
    
    
    
    
    
!------------!
!    损伤    !
!------------!
subroutine data_output_damage()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    real(8)::a
    integer i
   
    write( cStr , '(i5)' ) ii
    open( 10 , File = 'point.Damage' // AdjustL(Trim( cStr ) ) )
    write(10,'(a)')"Damage"
    write(10,'(a)')"part"
    write(10,'(a)')"         1"
    write(10,'(a)')"coordinates"  
	do i=1,totout(iload)
        ni=outnode(iload,i)
        a=damage(ni)    
        if(a.ge.0) then
            cFmt="(1x,e11.5)"
        else
            cFmt="(e12.5)"
        end if
        write( 10 , cFmt ) a
    end do 
    close(10)
end subroutine data_output_damage
!------------!
!    位移    !
!------------!
subroutine data_output_Displacement()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Displacement-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Displacement-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
            ni=outnode(iload,i)
            a=disp(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
end subroutine data_output_Displacement
!------------!
!    速度    !
!------------!
subroutine data_output_Velocity()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Velocity-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Velocity-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
                ni=outnode(iload,i)
            a=vel(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
end subroutine data_output_Velocity
 
subroutine data_output_Bforce()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Bforce-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Bforce-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
            ni=outnode(iload,i)
            a=bforce(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
end subroutine data_output_Bforce
    
subroutine data_output_Pforce()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Pforce-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Pforce-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
            ni=outnode(iload,i)
            a=Pforce(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
end subroutine data_output_Pforce

subroutine data_output_Cforce()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Cforce-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Cforce-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
            ni=outnode(iload,i)
            a=cforce(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
end subroutine data_output_Cforce   

subroutine data_output_Force()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Force-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Force-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
            ni=outnode(iload,i)
            a=Force(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
    end subroutine data_output_Force 
    
    subroutine data_output_Gforce()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    character  (len=20):: cstr
    character  (len=20):: cFmt
    character  (len=1):: xyz
    real(8)::a
    integer i
    do idime=1,ndime
        if (idime==1)then
             xyz='X'
        elseif(idime==2)then
             xyz='Y'
        elseif(idime==3)then
             xyz='Z'
        end if   
        write( cStr , '(i5)' ) ii
        open( 10 , File = 'point.Gforce-' //xyz// AdjustL(Trim( cStr ) ) )
        write(10,'(a)')"Gforce-"//xyz
        write(10,'(a)')"part"
        write(10,'(a)')"         1"
        write(10,'(a)')"coordinates" 
	    do i=1,totout(iload)
            ni=outnode(iload,i)
            a=gforce(ni,idime)     
            if(a.ge.0) then
                cFmt="(1x,e11.5)"
            else
                cFmt="(e12.5)"
            end if
            write( 10 , cFmt ) a
        end do 
        close(10)
    end do
end subroutine data_output_Gforce