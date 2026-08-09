# Use rgl null device so that rgl (a transitive dependency via fsbrain) does
# not attempt to load OpenGL/GLU system libraries on headless CI machines.
options(rgl.useNULL = TRUE)
