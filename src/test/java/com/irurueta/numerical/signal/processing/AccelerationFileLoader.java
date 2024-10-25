/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.irurueta.numerical.signal.processing;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

/**
 * Utility class to load .dat files containing accelerometer samples obtained on
 * a Nexus 6 device for testing purposes.
 */
class AccelerationFileLoader {

    static Data load(final File f) throws IOException {
        final var bufferSize = computeBufferSizes(f);
        final var stream = new DataInputStream(Files.newInputStream(f.toPath()));

        final var accelerationX = new float[bufferSize];
        final var accelerationY = new float[bufferSize];
        final var accelerationZ = new float[bufferSize];

        final var timestamp = new long[bufferSize];
        final var count = new long[bufferSize];

        for (var i = 0; i < bufferSize; i++) {
            count[i] = stream.readLong();
            timestamp[i] = stream.readLong();

            accelerationX[i] = stream.readFloat();
            accelerationY[i] = stream.readFloat();
            accelerationZ[i] = stream.readFloat();
        }

        stream.close();

        final var data = new Data();
        data.accelerationX = accelerationX;
        data.accelerationY = accelerationY;
        data.accelerationZ = accelerationZ;
        data.timestamp = timestamp;
        data.count = count;
        data.numSamples = bufferSize;

        return data;
    }

    private static int computeBufferSizes(final File f) {
        final var fileLength = f.length();

        final var floatSize = Float.SIZE / Byte.SIZE;
        final var longSize = Long.SIZE / Byte.SIZE;

        return (int) fileLength / (3 * floatSize + 2 * longSize);
    }

    public static class Data {
        float[] accelerationX;
        float[] accelerationY;
        float[] accelerationZ;
        long[] timestamp;
        long[] count;
        public int numSamples;
    }
}
