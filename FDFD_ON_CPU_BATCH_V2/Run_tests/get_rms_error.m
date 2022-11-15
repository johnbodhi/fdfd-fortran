function [rms_error] = get_rms_error(ydata, stt)
rms_error = (1/181)*sum(sqrt((ydata(1:181) - stt(1:181,2).').^2));
