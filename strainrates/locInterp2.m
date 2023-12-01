function [intXVel,intYVel] = locInterp2(rowCoord,colCoord,sqVx,sqVy)

ULxVel = sqVx(floor(rowCoord),floor(colCoord));
URxVel = sqVx(floor(rowCoord),ceil(colCoord));
LLxVel = sqVx(ceil(rowCoord),floor(colCoord));
LRxVel = sqVx(ceil(rowCoord),ceil(colCoord));

ULyVel = sqVy(floor(rowCoord),floor(colCoord));
URyVel = sqVy(floor(rowCoord),ceil(colCoord));
LLyVel = sqVy(ceil(rowCoord),floor(colCoord));
LRyVel = sqVy(ceil(rowCoord),ceil(colCoord));

topInterp = colCoord - floor(colCoord);
sideInterp = rowCoord - floor(rowCoord);

topxVel = ULxVel + (URxVel-ULxVel)*topInterp;
botxVel = LLxVel + (LRxVel-LLxVel)*topInterp;

topyVel = ULyVel + (URyVel-ULyVel)*topInterp;
botyVel = LLyVel + (LRyVel-LLyVel)*topInterp;

intXVel = topxVel + (botxVel-topxVel)*sideInterp;
intYVel = topyVel + (botyVel-topyVel)*sideInterp;

