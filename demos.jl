using HampelOutliers, Plots, LaTeXStrings

t = (1:50) / 10
x = [1:2:40; 5t + (@. 6cos(t + 0.5(t)^2)); fill(40,20)]
x[12] = -10
x[50:52] .= -12
x[79:82] .= [-5, 50, 55, 0]

plot(x, m=2, msw=0, label="Original", xlabel=L"k", ylabel=L"x_k")
m = Hampel.filter(x, 4, threshold=0)
scatter!(m, m=2, label="Median filter")
y = Hampel.filter(x, 4)
scatter!(y, m=2, label="Hampel filter")

savefig("demo_comparison.svg")
