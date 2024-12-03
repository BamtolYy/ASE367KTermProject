function guststate_new = dryden_gust_sample(gust_state,V,h,gustintensity,dt)
%gust state: [u_g,v_G,w_g] in m/s
%V: total airspeed in m/s
%h: altitude above ground in m
%gustintensity: gust_rms in m/s (1.5 m/s light, 3 m/s moderate, 4.5 m/s severe)
%dt: time step in s
whitenoise = randn([3,1]);
gustu = gust_state(1);
gustv = gust_state(2);
gustw = gust_state(3);
%rotated such that x vector is pointed directly upwards
hft = h*3.28084;
hft = max(hft,1);
if hft<1000
sigmaw = gustintensity/((0.177+0.000823*hft).^(0.4));
sigmav = sigmaw;
sigmau = gustintensity;
Lw = (hft/(0.177+0.000823*hft).^1.2)/3.28084;
Lv = Lw;
Lu = hft/3.28084;
elseif hft>2000
sigmaw = gustintensity;
sigmau = gustintensity;
sigmav = gustintensity;
Lu = 1750/3.28084;
Lv = Lu;
Lw = Lu;
else
sigmaw = gustintensity;
sigmau = gustintensity;
sigmav = gustintensity;
Lu = (1000 + 0.75*(hft-1000))/3.28084;
Lv = Lu;
Lw = Lu;
end
au = V/Lu;
Cu = exp(-au*dt);
av = V/Lv;
Cv = exp(-av*dt);
aw = V/Lw;
Cw = exp(-aw*dt);
%simplified form from MIL-F-8785C
gustu = Cu*gustu + sqrt((1-Cu.^2))*sigmau*whitenoise(1);
gustv = Cv*gustv + sqrt((1-Cv.^2))*sigmav*whitenoise(2);
gustw = Cw*gustw + sqrt((1-Cw.^2))*sigmaw*whitenoise(3);
guststate_new = [gustu;gustv;gustw];
end