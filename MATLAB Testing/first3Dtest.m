% Testing sim3D Unreal Representation for Thesis Progress

% Generate updatable world
world = sim3d.World('Output', @actorInputImpl);

% Create movable actor
box1 = sim3d.Actor(ActorName = 'Box1', Mobility = sim3d.utils.MobilityTypes.Movable);

% Make it a box with given dimensions
createShape(box1,'box',[0.25, 0.25, 0.25]);

% Give it a color
box1.Color = [0.59, 0.59, 0.59];

% Include it in the world
add(world, box1);

% Set viewer window POV
viewport = createViewport(world, Translation = [-2, 0, 0.5], Rotation = [0 -0.25 0]);

sampleTime = 0.01;
stopTime = 10;
run(world, sampleTime, stopTime);


% Get rid of it all
delete(world)

function actorInputImpl(world)

    % Sets the actor outputs (e.g. actor position to follow a path)
    % world.Actors.Box1.Translation(3) = world.Actors.Box1.Translation(3) + 0.01;
    % world.Actors.Box1.Translation(2) = world.Actors.Box1.Translation(2) + 0.02;
    world.Actors.Box1.Rotation(3) = world.Actors.Box1.Rotation(3) + 0.01;
    % world.Actors.Box1.Scale = world.Actors.Box1.Scale + [0.05 0 0.05];
    % world.Actors.Box1.Color = world.Actors.Box1.Color + [0 0 0.1];

end

