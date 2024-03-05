library('dorothea')
library(dplyr)

data = get(data('dorothea_mm', package="dorothea"))

# get target genes for tf
trans_targets_nfkb1 = data[data$tf == 'NFKB1', ]
trans_targets_nfkb2 = data[data$tf == 'NFKB2', ]
print(trans_targets_nfkb2)

# get tf for genes that rise
tf_cytokine = data[data$target == 'Ccl2', ]

Mode <- function(data) {
  ux <- unique(data)
  ux[which.max(tabulate(match(data, ux)))]
}
