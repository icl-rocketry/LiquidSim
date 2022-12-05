using Modia
using CoolProp
usePlotPackage("PyPlot")

Connector = Model( P = potential, h = potential, mdot = flow )

Part = Model( inlet = Connector, outlet = Connector, equations = :[
        P_inlet = inlet.P
        P_outlet = outlet.P
        h_inlet = inlet.h
        h_outlet = outlet.h
        0 = inlet.mdot + outlet.mdot
        mdot = inlet.mdot ] 
)

Restriction = Part | Model(C_d = 1.0, A=1.0, F_L = 1.0, substance = "N2O", equations = :[
    T = CoolProp.PropsSI("T","H",h_inlet,"P",P_inlet,substance)
    P_sat = CoolProp.PropsSI("T","H",h_inlet,"P",P_inlet,substance)
    rho = CoolProp.PropsSI("D","H",h_inlet,"P",P_inlet,substance)
    P_2 = max(P_outlet,P_inlet*(1-F_L^2)+P_sat*F_L^2)
    mdot = C_d * A * sqrt(2*rho*(P_inlet-P_2))
    h_inlet = h_outlet
    ]
)

Ambient = Model( P_ambient = 1.0, inlet = Connector, equations = :[
        inlet.P = P_ambient ] )

Tank = Model( P_internal = 1.0, h_internal = 1.0, substance = "N2O", outlet = Connector, equations = :[
    outlet.P = P_internal
    outlet.h = h_internal ] )

FeedSystemOxidizer = Model(
    substance = "N2O",
    tank = Tank | Map(P_internal=4e6, h_internal = CoolProp.PropsSI("H","P",4e6,"T",290,"N2O"), substance =:substance),
    valve = Restriction | Map(C_d = 0.7, A=pi*(1e-2)^3, F_L = 0.8, substance =:substance),
    injector = Restriction | Map(C_d = 0.7, A=pi*(3e-3)^2, F_L = 0.95, substance =:substance, P_inlet=Var(start=3.9e6), h_inlet=Var(start=3.9e5)),
    outlet = Connector,
    connect = :[
      (tank.outlet, valve.inlet)
      (valve.outlet, injector.inlet)
      (injector.outlet, outlet)
    ]
)

FeedSystemFuel = Model(
    substance = "METHANOL",
    tank = Tank | Map(P_internal=4e6, h_internal = CoolProp.PropsSI("H","P",4e6,"T",290,"METHANOL"), substance =:substance),
    valve = Restriction | Map(C_d = 0.7, A=pi*(1e-2)^3, F_L = 0.8, substance =:substance),
    injector = Restriction | Map(C_d = 0.7, A=pi*(3e-3)^2, F_L = 0.95, substance =:substance, P_inlet=Var(start=3.9e6), h_inlet=Var(start=3.9e5)),
    outlet = Connector,
    connect = :[
      (tank.outlet, valve.inlet)
      (valve.outlet, injector.inlet)
      (injector.outlet, outlet)
    ]
)

UpperEngine = Model(
    fuel = "METHANOL",
    oxidizer = "N2O",
    #T_tank_oxidizer = 290,
    #P_tank_oxidizer = 4e6;
    #T_tank_fuel = 290,
    #P_tank_fuel = 4e6;
    feed_system_fuel = FeedSystemFuel | Map(substance =:fuel),
    feed_system_oxidizer = FeedSystemOxidizer | Map(substance =:oxidizer),
    ambient_fuel = Ambient | Map(P_ambient=2e6),
    ambient_oxidizer = Ambient | Map(P_ambient=2e6),
    connect = :[
      (feed_system_fuel.outlet, ambient_fuel.inlet)
      (feed_system_oxidizer.outlet, ambient_oxidizer.inlet)
    ]
)

upper_engine = @instantiateModel(UpperEngine)
simulate!(upper_engine, stopTime=10, merge=Map(T=0.5, x=0.8))