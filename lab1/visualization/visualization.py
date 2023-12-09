import matplotlib.pyplot as plt
from celluloid import Camera

pointsNumber = 20

def readFromfile(filePath, pointsNumber):
    points = [[] for _ in range(pointsNumber)]
    with open(filePath, 'r') as file:
        for line in file:
            lineParts = line.strip().split()
            if len(lineParts) > 2:
                pointNumber = int(lineParts[1])
                x = lineParts[3]
                y = lineParts[4]
                points[pointNumber].append({"x": float(x), "y": float(y)})
    return points


points = readFromfile("output.txt", pointsNumber)

fig, ax = plt.subplots()
camera = Camera(fig)

for i in range(len(points[0])):
    for point in points:
        plt.scatter(point[i]["x"], point[i]["y"])
    camera.snap()
    
animation = camera.animate()
animation.save('point_movement.gif')

plt.show()