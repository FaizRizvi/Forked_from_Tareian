#Colors
flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
colorful = sns.xkcd_palette(colors)
grade = sns.light_palette("green")
diverging = sns.diverging_palette(10, 220, sep=80, n=7)
green = sns.cubehelix_palette(8, start=.5, rot=-.75)
cube = sns.color_palette("cubehelix", 8)
pink = sns.cubehelix_palette(8)
teal = sns.color_palette("GnBu_d")
husl = sns.color_palette("husl", 8)

df = df.sort_values("Peaks")




u = sns.barplot(x="Sample", y="Peaks", data=df, palette=green)
for item in u.get_xticklabels():
    item.set_rotation(60)
u.set_title("Peaks Called Per ENCODE Library")

