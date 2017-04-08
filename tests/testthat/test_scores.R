messageU("scores")

data("combined_mode-differential")
info0 <- combined_model$info
bins0 <- combined_model$bins[1:9]
bins0$peakScores <- matrix(1:72, ncol=nrow(info0))
colnames(bins0$peakScores) <- info0$ID

### Test first case (normal) ###
bins <- bins0
info <- info0
mat <- bins$peakScores
differential.score <- differentialScoreSum(mat = mat, info = info)
expect_that(differential.score, equals(rep(90, 9)))

### Same mark, different conditions ###
mask <- info0$mark == 'H3K27me3'
bins <- bins0
info <- info0[mask,]
bins$peakScores <- bins$peakScores[,mask]
mat <- bins$peakScores
differential.score <- differentialScoreSum(mat = mat, info = info)
expect_that(differential.score, equals(rep(45, 9)))

### Same condition, different marks ###
mask <- info0$condition == 'BN'
bins <- bins0
info <- info0[mask,]
bins$peakScores <- bins$peakScores[,mask]
mat <- bins$peakScores
differential.score <- differentialScoreSum(mat = mat, info = info)
expect_that(differential.score, equals(rep(18, 9)))

### Test first case (normal) with only 1 bin ###
bins <- bins0[1]
info <- info0
mat <- bins$peakScores
differential.score <- differentialScoreSum(mat = mat, info = info)
# expect_that(differential.score, equals(54))
expect_that(differential.score, equals(90))

### Same mark, different conditions with only 1 bin ###
mask <- info0$mark == 'H3K27me3'
bins <- bins0[1]
info <- info0[mask,]
bins$peakScores <- bins$peakScores[,mask,drop=FALSE]
mat <- bins$peakScores
differential.score <- differentialScoreSum(mat = mat, info = info)
expect_that(differential.score, equals(45))

### Same condition, different marks with only 1 bin ###
mask <- info0$condition == 'BN'
bins <- bins0[1]
info <- info0[mask,]
bins$peakScores <- bins$peakScores[,mask,drop=FALSE]
mat <- bins$peakScores
differential.score <- differentialScoreSum(mat = mat, info = info)
expect_that(differential.score, equals(18))
