round_smart <- function(x,place = 1,direction="nearest"){

  if(direction%in%c("up","UP","Up","uP")){direction = "up"}
  if(direction%in%c("down","Down","DOWN")){direction = "down"}

  res = switch(direction,
               "nearest" = round(x/place)*place,
               "up" = ceiling(x/place)*place,
               "down" = floor(x/place)*place)

  return(res)
}
