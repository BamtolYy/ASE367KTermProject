function [state,aux] = vehicle_dynamics(state,params,dt)
    x = state(1);
    z = state(2);
    vx = state(3);
    vz = state(4);
    theta = state(5);
    q = state(6);

    gust_state = state(7:9);

    mass = params(1);
    Ixx = params(2);
    A = params(3);
    cna = params(4);
    cnq = params(5);
    cma = params(6);
    cmq = params(7);
    cd = params(8);
    gustintensity = params(9);
    gimbalcg = params(10);
    
    R = 6371000;
    mu = 3.986e14;
    T_max = 15500;

    fpa = atan2(vx,vz);

    r = sqrt((z+R).^2+(x).^2);
    h = r-R;
    psi = atan2(x,z+R);

    rho = 1.225*exp(-h/10400);

    g = mu/r.^2;

    gust_state = dryden_gust_sample(gust_state,sqrt(vx.^2+vz.^2),h,gustintensity,dt);
    gust_x = gust_state(1);
    gust_w = gust_state(2);

    u = vz*cos(theta) + vx*sin(theta)-gust_x;
    w = vx*cos(theta) - vz*sin(theta)-gust_w;

    alpha = atan2(-w,u);

    Q = 0.5*rho*(u.^2+w.^2);

    N = Q*A*(cna*alpha+cnq*q);
    D = Q*A*cd;

    [delta_gimbal,delta_T,control_theta] = rocket_controller(state,params);

    T = T_max*delta_T;

    Fx = T*cos(delta_gimbal) - D;
    Fz = T*sin(delta_gimbal) + N;

    ax = (Fz*cos(theta)+Fx*sin(theta))/mass - g*sin(psi);
    az = (Fx*cos(theta)-Fz*sin(theta))/mass - g*cos(psi);

    M = Q*A*(cma*alpha)/Ixx + (exp(Q*A*cmq/Ixx*dt)-1)*q/dt - gimbalcg*sin(delta_gimbal)*T/Ixx;

    M1 = Q*A*(cma*alpha)/Ixx;
    M2 = (exp(Q*A*cmq/Ixx*dt)-1)*q/dt;

    dx = [vx+0.5*ax*dt;vz+0.5*az*dt; ax; az; q; M; 0; 0; 0];

    state = state+dx*dt;
    state(7:9) = gust_state;
    aux = [alpha,rho,Q,ax,az,psi,fpa,delta_gimbal,control_theta,delta_T];
    
end