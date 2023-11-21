# """ Workshop 9: trait hypervolumes
#     BSC 6926 B53
#     date: 11/20/2023"""


# This workshop discusses working with community trait data by introducing the `hypervolume` package by going over how to prepare the data and how to build hypervolumes in multiple ways.\
# 
# [abundance data](https://github.com/BSC6926-B52-Fall2023/workshopScripts/blob/main/data/abundHermine.csv)
# 
# [trait data](https://github.com/BSC6926-B52-Fall2023/workshopScripts/blob/main/data/fishTraits.csv)

# ## Merge/Join
# If two data frames contain different columns of data, then they can be merged together with the family of join functions.
# 
# +`left_join()` = uses left df as template and joins all matching columns from right df 
# +`right_join()` = uses right df as template and joins all matching columns from left df
# +`inner_join()` = only matches columns contained in both dfs
# +`full_join()` = combines all rows in both dfs

library(tidyverse)

left = tibble(name = c('a', 'b', 'c'),
              n = c(1, 6, 7), 
              bio = c(100, 43, 57))

right = tibble(name = c('a', 'b', 'd', 'e'),
               cals = c(500, 450, 570, 600))

left_join(left, right, by = 'name')

right_join(left, right, by = 'name')

inner_join(left, right, by = 'name')

full_join(left, right, by = 'name')

# multiple matches
fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))

col = tibble(species = c('Salmon', 'Cod'),
             coast = c('West', 'East'))

left_join(fish, col, by = 'species')


# ## scaling data
# Because hypervolumes can be generated with any continuous data as an axes, many of the times the units are not combatible. Blonder et al. [2014](https://doi-org.ezproxy.fiu.edu/10.1111/geb.12146) & [2018](https://doi-org.ezproxy.fiu.edu/10.1111/2041-210X.12865) to convert all of the axes into the same units. This can be done by taking the z-score of the values to convert units into standard deviations. Z-scoring data can be done with the formula:
#   $$ z = \frac{x_{i}-\overline{x}}{sd} $$ Where $x_{i}$ is a value, $\overline{x}$ is the mean, and $sd$ is the standard deviation. By z-scoring each axis, 0 is the mean of that axis, a value of 1 means that the value is 1 standard deviation above the global mean of that axis, and a value of -1 is 1 standard deviation below the global mean of the axis. In R this can be done manually or with the `scale()` function. 

fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))

#
fish = fish |> 
  mutate(zcatch1 = (catch - mean(catch))/sd(catch), # manual
         zcatch2 = scale(catch)) # with scale

fish 

# center = mean, scale = sd
fish$zcatch2


# ## nesting data
# One benefit of `tibbles` is that they can contain list columns. This means that we can make columns of `tibbles` that are nested within a dataset. Nesting creates a list-column of data frames; unnesting flattens it back out into regular columns. Nesting is a implicitly summarising operation: you get one row for each group defined by the non-nested columns. This is useful in conjunction with other summaries that work with whole datasets, most notably models. This can be done with the `nest()` and then flattened with `unnest()`

fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))

# using group_by
fish_nest = fish |> 
  group_by(species) |> 
  nest()

fish_nest
fish_nest$data

# using .by in nest
# column name becomes data unless you change .key
fish_nest2 = fish |> 
  nest(.by = year, .key = 'df')

fish_nest2
fish_nest2$df

## map
### `purr`
# The newest and new standard package with `tidyverse` is `purr` with its set of `map()` functions. Some similarity to `plyr` (and base) and `dplyr` functions but with more consistent names and arguments. Notice that map function can have some specification for the type of output.
# + `map()` makes a list.
# + `map_lgl()` makes a logical vector.
# + `map_int()` makes an integer vector.
# + `map_dbl()` makes a double vector.
# + `map_chr()` makes a character vector.

df = iris  |> 
  select(-Species)
#summary statistics
map_dbl(df, mean)

# using map with mutate and nest
d = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
           year = rep(c(1999,2005,2020), each = 2),
           catch = c(50, 60, 40, 50, 60, 100)) |> 
  nest(.by = species) |> 
  mutate(correlation = map(data, \(data) cor.test(data$year, data$catch)))

d
d$correlation

# ## Hypervolumes 
# Hypervolumes are a multidimensional tool that is based on Hutchinson's *n*-dimensional niche concept and we can build them with the `hypervolume` package. 
# 
# ### Preparing the data
# Typically we have a dataset that has the abundance of the community and another that has trait data. Therefore we need to combine them. Also, we can't have `NAs`, so we have to filter any missing data one we combine the datasets. 

# abundance data
ab = read_csv('data/abundHermine.csv')

# trait data
tr = read_csv('data/fishTraits.csv')

# combine
df = left_join(ab, tr, by = 'species') |> 
  drop_na()

# Once combined we can select only the columns we want in the hypervolume, z-score and nest the data.

df = df |> 
  select(period, trophic_level, temp_preference, generation_time) |> 
  mutate(across(trophic_level:generation_time, scale)) |> 
  group_by(period) |> 
  nest()

df

### Building hypervolumes
# With a nested dataset of our columns that we want to build hypervolumes for we can use `mutate()` and `map()` to generate the hypervolume. 

library(hypervolume)

df = df |> 
  mutate(hv = map(data, ~hypervolume_gaussian(.x, name = period,
                                              samples.per.point = 1000,
                                              kde.bandwidth = estimate_bandwidth(.x), 
                                              sd.count = 3, 
                                              quantile.requested = 0.95, 
                                              quantile.requested.type = "probability", 
                                              chunk.size = 1000, 
                                              verbose = F)))

df
df$hv
# ### plotting hypervolumes 
# We can plot multiple hypervolumes by joining them together 

hvj = hypervolume_join(df$hv[[1]], df$hv[[2]])

plot(hvj, pairplot = T, colors=c('goldenrod','blue'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

# ### hypervolume metrics
# The geometry of hypervolumes are useful when characterizing and comparing hypervolumes. Hypervolume size represents the variation of the data, centroid distance compares the euclidian distance between two hypervolume centroids (mean conditions), and overlap measures the simularity of hypervolumes. 

# size 
df = df |> 
  mutate(hv_size = map_dbl(hv, \(hv) get_volume(hv)))

df

# centroid distance 
hypervolume_distance(df$hv[[1]], df$hv[[2]], type = 'centroid', check.memory=F)

# overlap 
hvset = hypervolume_set(df$hv[[1]], df$hv[[2]], check.memory = F)

hypervolume_overlap_statistics(hvset)

# ### Weight hypervolume input
# The above hypervolume is just based on the traits using presence of species, but we can weight the points to shape the hypervolume based on abundance 

#prep data
df_w = left_join(ab, tr, by = 'species') |> 
  drop_na() |> 
  select(period, abund, trophic_level, temp_preference, generation_time) |> 
  mutate(across(trophic_level:generation_time, scale)) |> 
  group_by(period) |> 
  nest(weight = abund, data = trophic_level:generation_time) 
df_w

# make hypervolumes
df_w = df_w |> 
  mutate(hv = map2(data,weight, ~hypervolume_gaussian(.x, name = paste(period,'weighted',sep = '_'),
                                                      weight = .y$abund,
                                                      samples.per.point = 1000,
                                                      kde.bandwidth = estimate_bandwidth(.x), 
                                                      sd.count = 3, 
                                                      quantile.requested = 0.95, 
                                                      quantile.requested.type = "probability", 
                                                      chunk.size = 1000, 
                                                      verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)))

df_w


hvj_w = hypervolume_join(df_w$hv[[1]], df_w$hv[[2]])

plot(hvj_w, pairplot = T, colors=c('goldenrod','blue'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

# centroid distance 
hypervolume_distance(df_w$hv[[1]], df_w$hv[[2]], type = 'centroid', check.memory=F)

# overlap 
hvset_w = hypervolume_set(df_w$hv[[1]], df_w$hv[[2]], check.memory = F)

hypervolume_overlap_statistics(hvset_w)

# ### From mean and sd
# Sometimes, we do not have enough points to meet assumptions to make a hypervolume. Therefore, we can simulate random points based on mean and sd of our axes. We can then simulate the information needed and make our hypervolumes.

# mean and sd
df_m = left_join(ab, tr, by = 'species') |> 
  drop_na() |> 
  pivot_longer(trophic_level:generation_time, names_to = 'trait', values_to = 'value') |> 
  group_by(period,trait) |> 
  summarize(mean = mean(value),
            sd = sd(value))

df_m

#generate points 
# number of points 
n = 50 

df_tot = df_m |> slice(rep(1:n(), each=n))|>
  mutate(point = map2_dbl(mean,sd, \(mean,sd) rnorm(1,mean =mean,sd =sd))) |> 
  group_by(period, trait) |> 
  mutate(num = row_number()) |>
  select(-mean, -sd)|>
  pivot_wider(names_from = trait, values_from = point)|> 
  select(-num) |> 
  mutate(across(generation_time:trophic_level,scale)) |> 
  group_by(period) |> 
  nest()


# generate hypervolumes
df_tot = df_tot |> 
  mutate(hv = map(data, ~hypervolume_gaussian(.x, name = period,
                                              samples.per.point = 1000,
                                              kde.bandwidth = estimate_bandwidth(.x), 
                                              sd.count = 3, 
                                              quantile.requested = 0.95, 
                                              quantile.requested.type = "probability", 
                                              chunk.size = 1000, 
                                              verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)))


#plot
hvj_tot = hypervolume_join(df_tot$hv[[1]], df_tot$hv[[2]])

plot(hvj_tot, pairplot = T, colors=c('goldenrod','blue'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

# centroid distance 
hypervolume_distance(df_tot$hv[[1]], df_tot$hv[[2]], type = 'centroid', check.memory=F)

# overlap 
hvset_tot = hypervolume_set(df_tot$hv[[1]], df_tot$hv[[2]], check.memory = F)

hypervolume_overlap_statistics(hvset_tot)



