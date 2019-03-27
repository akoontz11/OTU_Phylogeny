# %%% Delta transformations example %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 <- sim.bdtree(b=0.5, d=0.2, stop="time", t=5)
# Removing extinct species
p1 <- drop.extinct(p1)
plot(p1, main="Test Tree")
class(p1)

# Examples of delta transformations:
recent.tree <- rescale(p1, "delta", 3)
plot(recent.tree, main="Trait evolution sped up over time (delta = 3)")

ancient.tree <- rescale(p1, "delta", 0.33)
plot(ancient.tree, main="Trait evolution slowed down over time (delta = 0.3)")
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%