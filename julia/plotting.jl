time = collect(0:19)
c = ones(19,1)*2.2
uu = [2.2;2.16353;2.008;2.0027;2.0009;2.00032;2.00013;2.00006;2.00004;2.00004;2.00004;2.00003;2.00003;2.00003;2.00003;2.00003;2.00003;2.00003;2.00003]
legend = ["constraints" "control action"]
time = [time time]
p = [c uu]
plot(time,p,
  line=[:steppost :steppost],
  lab=map(string,legend),
  linestyle=[:dash :solid],
  ylims = (1.8,2.4),
  xlims = (0,18),
  color=[RGB(0.87,0.49,0) RGB(0,0.45,0.74)],
  xlabel = "t[s]",
  ylabel = "u",
  linewidth = 3,
  xticks = 0:2:18,
  yticks = 1.8:0.1:2.4,
  title = "Control Action")

time = collect(0:19)
r = ones(19,1)*2
yy = [0;-0.52;0.7815;1.6037;1.871;1.958;1.986;1.9955;1.9985;1.9995;1.9998;1.99995;1.99998;1.99999;2;2;2;2;2]
legend = ["reference" "output"]
time = [time time]
p = [r yy]
plot(time,p,
  line=[:steppost :steppost],
  lab=map(string,legend),
  linestyle=[:dash :solid],
  ylims = (-1,2.5),
  xlims = (0,18),
  color=[RGB(0.75,0.75,0) RGB(0.85,0.33,0.1)],
  xlabel = "t[s]",
  ylabel = "y",
  linewidth = 3,
  xticks = 0:2:18,
  title = "MPC with Integral Control Action")
