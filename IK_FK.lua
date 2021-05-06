--------------------------------------------------
--
--  Author: Daniel Slovich
--  Date:   5/5/21
--  Info:   FK/IK Calculations
--          using Matrix Method
--          ECE 470 Humanoid Robotics
--          
--
--
--------------------------------------------------
matrix = require 'matrix'

--[[--Psuedo Inverse Functions--]]---------------------------------------------------------------
function pinv(a)
    local x = #a
    local y = #a[1]
    if x<y then
        b = pinv_v(a)
    elseif x>y then
        b = pinv_h(a)
    else
        b = matrix.invert(a)
    end
    return b
end

function pinv_h(a)
    at  = matrix.transpose(a)
    ata = matrix.mul(at,a)
    inv_ata = matrix.invert(ata)
    b = matrix.mul(inv_ata,at)
    return b
end

function pinv_v(a)
    at  = matrix.transpose(a)
    aat = matrix.mul(a, at)
    inv_aat = matrix.invert(aat)
    b = matrix.mul(at,inv_aat)
    return b
end

--[[--Jacobian Solver--]]------------------------------------------------------------------------
function jacobian(theta,dtheta)
    
    -- Pulling EEF exact location using simulation calls to avoid error propogation in FK.
    pt = sim.getObjectPosition(sim.getObjectHandle("L5"), -1, p)
    x0 = pt[1]
    y0 = pt[2]
    z0 = pt[3]
    
    dthetap ={dtheta,0,0,0,0}
    pt = get_fk(dthetap)
    pt1= get_fk1(dthetap)
    xdt1 = pt1[1]
    ydt1 = pt[2]
    zdt1 = pt1[2]
    
    dthetap ={0,dtheta,0,0,0}
    pt = get_fk(dthetap)
    pt1= get_fk1(dthetap)
    xdt2 = pt1[1]
    ydt2 = pt[2]
    zdt2 = pt1[2]
    
    dthetap ={0,0,dtheta,0,0}
    pt = get_fk1(dthetap)
    pt1= get_fk1(dthetap)
    xdt3 = pt1[1]
    ydt3 = pt[2]
    zdt3 = pt1[2]
    
    dthetap ={0,0,0,dtheta,0}
    pt = get_fk1(dthetap)
    pt1= get_fk1(dthetap)
    xdt4 = pt1[1]
    ydt4 = pt[2]
    zdt4 = pt1[2]
    
    dthetap ={0,0,0,0,dtheta}
    pt = get_fk1(dthetap)
    pt1= get_fk1(dthetap)
    xdte = pt1[1]
    ydte = pt[2]
    zdte = pt1[2]
    
    dxdt1 = (xdt1-x0)/dtheta
    dxdt2 = (xdt2-x0)/dtheta
    dxdt3 = (xdt3-x0)/dtheta
    dxdt4 = (xdt4-x0)/dtheta
    dxdte = (xdte-x0)/dtheta

    dydt1 = (ydt1-y0)/dtheta
    dydt2 = (ydt2-y0)/dtheta
    dydt3 = (ydt3-y0)/dtheta
    dydt4 = (ydt4-y0)/dtheta
    dydte = (ydte-y0)/dtheta

    dzdt1 = (zdt1-z0)/dtheta
    dzdt2 = (zdt2-z0)/dtheta
    dzdt3 = (zdt3-z0)/dtheta
    dzdt4 = (zdt4-z0)/dtheta
    dzdte = (zdte-z0)/dtheta
    
    J = matrix{ {dxdt1, dxdt2, dxdt3, dxdt4, dxdte},
                {dydt1, dydt2, dydt3, dydt4, dydte},
                {dzdt1, dzdt2, dzdt3, dzdt4, dzdte}}
                
    return J
end


--[[--Forward Kinematics Using Matrix Method--]]------------------------------------------------
function get_fk(th)                                     --  th = joint theta
    len  = {0.3,0.3,0.3,0.3,0.06}                       -- link length
    for i=1,5,1 do
        R[i] = {{math.cos(th[i]), -math.sin(th[i]), 0}, -- Find rotational mtx's
                {math.sin(th[i]),  math.cos(th[i]), 0},
                {0,               0,                1}}
                
        P[i] = {{1, 0,    0},                           -- Find positional mtx's
                {0, 1, len[i]},
                {0, 0,    1}}
        
        T[i] = matrix.mul(R[i], P[i])                   -- Translational mtx's for each link
        
        -- This should hopfully minimize error propogation in matrix.mul
        for j=1,3,1 do
            for k =1,3,1 do
              if (T[i][j][k]<0.0001 and T[i][j][k]>-0.0001) then
                  T[i][j][k] = 0
              end
            end
         end 
 
    end

    T_01 = T[1]                                          -- Overall Translational mtx's
    T_02 = matrix.mul(T_01,T[2])
    T_03 = matrix.mul(T_02,T[3])
    T_04 = matrix.mul(T_03,T[4])
    T_05 = matrix.mul(T_04,T[5])
    fk = {T_05[1][3],T_05[2][3],T_05[1][2]}
    return(fk)
end

--[[--Forward Kinematics Using Analytic Method--]]---------------------------------------------
function get_fk1(th)                                    -- th = joint theta
    x ={0,0,0,0,0}                                      -- init x
    y ={0,0,0,0,0}                                      -- init y
    z ={0.0,0,0,0,0}                                    -- init z
    len  = {0.3,0.3,0.3,0.3,0.06} 
    Theta = th[1]                                       -- init theta sum
    
    for i=2,5,1 do                                      -- find (x,z) at each joint
        Theta = Theta + (-th[i])
        x[i]=x[i-1] + len[i]*math.sin(Theta)
        z[i]=z[i-1] + len[i]*math.cos(Theta)
    end
                                      
    return {x[5], 0, z[5], Theta}
end

--[[--Initialization Function--]]-----------------------------------------------------------
function sysCall_init()                                 -- SYS_INIT Function
    j1=sim.getObjectHandle("J1")                        -- Joint Object Handles
    j2=sim.getObjectHandle("J2")
    j3=sim.getObjectHandle("J3")
    j4=sim.getObjectHandle("J4")
    j5=sim.getObjectHandle("J5")
    
    scale = 0.001;                                      -- dynamic scaling smooths movement
    move = 0.05                                         -- init val needed to push block
    t  = {0.0,0,0,0,0}                                  -- theta
    j  = {j1,j2,j3,j4,j5}                               -- joints
    jp = {}                                             -- joint pos
    R  = {}                                             -- Rotational mtx's MTX
    P  = {}                                             -- Position mtx's MTX
    T  = {}                                             -- Translational mtx's MTX
    b  = 0.05                                           -- base height
end

--[[--Actuation Function--]]--------------------------------------------------------------
function sysCall_actuation()                            -- SYS_ACTUATE Function
    move = move-0.00001
    dtheta  = 0.001
    -- Gather exact goal and eef locations using sim calls to avoid error propogation in FK
    goal    = sim.getObjectPosition(sim.getObjectHandle("Midterm"), -1, p)
    current = sim.getObjectPosition(sim.getObjectHandle("L5"), -1, p)
    
    if scale >= 0.9005 then
        scale = 0.9005
    else
        scale = scale + 0.0005
    end
    
    for i=1,5,1 do                                      -- Sets joint pos to thetas
        sim.setJointTargetPosition(j[i], t[i])
    end
    
    for i=1,5,1 do
        jp[i] = sim.getJointPosition(j[i])              -- Gets joint pos thetas
    end
    
    if(goal[3]<=0.4 and goal[3]>=-0.4) then             -- Halts arm once block is pushed off
        t = t
    else
    
        error = matrix{ move+(goal[1]-current[1]),      -- error to goal/move pushes block
                        0.00+(goal[2]-current[2]),
                        0.00+(goal[3]-current[3])}
            
        J      = jacobian(jp, dtheta)                   -- Using Inv Jacobian to eval dthetas
        Jt     = matrix.transpose(J)
        invJ   = pinv(J)
        delta  = matrix.mul(invJ, error)
        
        t[1] = scale*0.4000*(jp[1]+delta[1][1])         -- New thetas with custom joint scaling
        t[2] = scale*1.0600*(jp[2]+delta[2][1])         -- Joint scaling makes for a smoother
        t[3] = scale*1.1100*(jp[3]+delta[3][1])         -- and more accurate movement to goal
        t[4] = scale*1.1095*(jp[4]+delta[4][1])
        t[5] = scale*1.0000*(jp[5]+delta[5][1])
    end
end

function sysCall_sensing()
    -- put your sensing code here
end

function sysCall_cleanup()
    -- do some clean-up here
end
