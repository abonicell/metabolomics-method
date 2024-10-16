hp_plot <- Hilic_Pos$data
hp_plot$class <- Hilic_Pos$sample_meta$class

hp_plot %>% 
  select(class, "alpha-Aminodiphenylacetic acid":"Xanthine") %>% 
  gather(Parameter, Value, -class) %>% 
  ggplot(aes(class, Value, fill = class)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_light() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          size = rel(0.8),
          angle = 45,
          hjust = 1,
          vjust = 1.1
        )
  ) +
  scale_fill_manual(values = c("#386cb0","#ef3b2c","#7fc97f","black"))


hn_plot <- Hilic_Neg$data
hn_plot$class <- Hilic_Neg$sample_meta$class

hn_plot %>% 
  select(class, "2-Ethyl-2-hydroxybutyric acid":"Undecanedioic acid" ) %>% 
  gather(Parameter, Value, -class) %>% 
  ggplot(aes(class, Value, fill = class)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_light() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          size = rel(0.8),
          angle = 45,
          hjust = 1,
          vjust = 1.1
        )
  ) +
  scale_fill_manual(values = c("#386cb0","#ef3b2c","#7fc97f","black"))


rpp_plot <- C18_Pos$data
rpp_plot$class <- C18_Pos$sample_meta$class

rpp_plot %>% 
  select(class, "alpha-Ergocryptine":"Val-Val-Arg") %>% 
  gather(Parameter, Value, -class) %>% 
  ggplot(aes(class, Value, fill = class)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_light() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          size = rel(0.8),
          angle = 45,
          hjust = 1,
          vjust = 1.1
        )
  ) +
  scale_fill_manual(values = c("#386cb0","#ef3b2c","#7fc97f","black"))


rpn_plot <- C18_Neg$data
rpn_plot$class <- C18_Neg$sample_meta$class

rpn_plot %>% 
  select(class, "(+/-)-5-Hydroxy-6E,8Z,11Z,14Z,17Z-eicosapentaenoic acid":"trans-2-Hydroxycinnamic acid") %>% 
  gather(Parameter, Value, -class) %>% 
  ggplot(aes(class, Value, fill = class)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_light() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          size = rel(0.8),
          angle = 45,
          hjust = 1,
          vjust = 1.1
        )
  ) +
  scale_fill_manual(values = c("#386cb0","#ef3b2c","#7fc97f","black"))
