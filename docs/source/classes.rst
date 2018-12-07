Classes within Pyplanes
=======================


FEM Package
-----------

.. uml::

  package fem.media <<Frame>> {
    abstract fem.media.FEMMedium {
      {abstract} get_physical_multipliers()
    }
    abstract fem.media.PEM

    fem.media.FEMMedium <|-- fem.media.Fluid
    fem.media.FEMMedium <|-- fem.media.Elastic
    fem.media.FEMMedium <|-- fem.media.PEM

    fem.media.PEM <|-- fem.media.PEM1998
    fem.media.PEM <|-- fem.media.PEM2001
  }

  package "media" <<Frame>> {
    abstract media.Medium {
      MODEL
      TYPE

      from_dict()
      {abstract} compute_missing()
      {abstract} update_frequency()
    }

    media.Medium <|-- media.Fluid
    media.Medium <|-- media.Elastic
    media.Medium <|-- media.EqFluidJCA
    media.EqFluidJCA <|-- media.PEM
    media.EqFluidJCA <|-- media.Limp
  }

  fem.media.Fluid <|.. media.Fluid
  fem.media.Fluid <|.. media.Limp







